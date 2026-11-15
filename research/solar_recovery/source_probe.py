from pathlib import Path
import json, hashlib, datetime, requests
O=Path("solar-evidence");O.mkdir(exist_ok=True)
S=requests.Session()
S.headers["User-Agent"]="Solar-Intelligence-academic-source-audit/0.9"
SER="hmi.sharp_cea_720s";T="2017.09.04_00:00:00_TAI"
FIELDS="HARPNUM,NOAA_AR,NOAA_NUM,NOAA_ARS,T_REC,T_OBS,QUALITY,CRPIX1,CRPIX2,CRVAL1,CRVAL2,CDELT1,CDELT2,CTYPE1,CTYPE2"
SEGS=["continuum","Br","Bp","Bt","bitmap"]
RESULT={"status":"REVIEW_REQUIRED_NOT_A_TRAINING_SET","harps":[],"noaa":[]}
def sha(b):return hashlib.sha256(b).hexdigest()
def out(filename,obj):
 p=O/filename;p.parent.mkdir(parents=True,exist_ok=True)
 p.write_text(json.dumps(obj,indent=2,default=str))
def get(url,params=None,cap=14000000):
 with S.get(url,params=params,timeout=(15,60),stream=True) as r:
  r.raise_for_status()
  body=bytearray()
  for piece in r.iter_content(131072):
   body+=piece
   if len(body)>cap:raise ValueError("FILE_EXCEEDS_LIMIT")
  return bytes(body),r.url
def query(ds,keys,segs=None):
 attempts=[]
 for host in ("https://jsoc.stanford.edu","https://jsoc1.stanford.edu"):
  try:
   q={"op":"rs_list","ds":ds,"key":keys}
   if segs:q["seg"]=segs
   data,url=get(host+"/cgi-bin/ajax/jsoc_info",q,3000000)
   doc=json.loads(data)
   if int(doc.get("status",-1))!=0:raise ValueError("DRMS_STATUS_"+str(doc.get("status")))
   return {"ok":True,"url":url,"sha256":sha(data),"raw_bytes":len(data),"data":doc,"failed_hosts":attempts}
  except Exception as e:attempts.append({"host":host,"error":type(e).__name__+":"+str(e)[:200]})
 return {"ok":False,"attempts":attempts}
def rows(d,kind):
 n=int(d.get("count",0))
 if n>1000:raise ValueError("TOO_MANY_ROWS")
 a=[{} for _ in range(n)]
 for col in d.get(kind,[]):
  if len(col.get("values",[]))!=n:raise ValueError("BAD_COL_"+kind)
  for i,v in enumerate(col["values"]):a[i][col["name"]]=v
 return a
def download(harp,seg,path):
 result={"segment":seg,"state":"SOURCE_UNAVAILABLE"}
 if not isinstance(path,str) or not path.startswith("/SUM") or ".." in path:
  return {"segment":seg,"state":"UNSAFE_OR_MISSING_SEGMENT_PATH","path":str(path)}
 for host in ("https://jsoc1.stanford.edu","https://jsoc.stanford.edu"):
  try:
   raw,url=get(host+path)
   if not raw.startswith(b"SIMPLE  =") or len(raw)%2880:
    raise ValueError("FITS_STRUCTURE_INVALID")
   d=O/("fits_"+str(harp));d.mkdir(exist_ok=True);(d/(seg+".fits")).write_bytes(raw)
   result.update(state="SOURCE_FITS_ACQUIRED",bytes=len(raw),sha256=sha(raw),url=url)
   try:
    from astropy.io import fits
    import numpy as np
    with fits.open(d/(seg+".fits"),checksum=True,memmap=False) as hdus:
     arrays=[np.asarray(x.data) for x in hdus if getattr(x.data,"ndim",0)==2]
     if len(arrays)==1:
      a=arrays[0];valid=np.isfinite(a)
      result["shape"]=list(a.shape);result["finite_fraction"]=float(valid.mean())
      result["range"]=[float(a[valid].min()),float(a[valid].max())] if valid.any() else None
     result["checksums"]=[x.verify_checksum() if "CHECKSUM" in x.header else None for x in hdus]
   except Exception as e:result["pixel_error"]=str(e)[:180]
   return result
  except Exception as e:result.setdefault("attempts",[]).append(type(e).__name__+":"+str(e)[:160])
 return result
reverse=query(f"{SER}[][{T}]","HARPNUM,NOAA_AR,NOAA_ARS,NOAA_NUM,T_REC,QUALITY")
out("reverse_cea.json",reverse)
for harp,noaa in ((7115,12673),(7122,12677)):
 r={"harp":harp,"candidate_noaa":noaa,"accepted":False}
 q=query(f"{SER}[{harp}][{T}]",FIELDS,",".join(SEGS))
 out("query_"+str(harp)+".json",q)
 r["source_metadata_ok"]=q["ok"];r["reverse_metadata_ok"]=reverse["ok"]
 if q["ok"]:
  try:
   rr=rows(q["data"],"keywords");sg=rows(q["data"],"segments")
   r["exact_rows"]=rr
   if reverse["ok"]:r["reverse_rows"]=rows(reverse["data"],"keywords")
   if len(rr)==len(sg)==1:r["files"]=[download(harp,s,sg[0].get(s)) for s in SEGS]
   else:r["error"]="NOT_SINGLE_CEA_ROW"
  except Exception as e:r["error"]=type(e).__name__+":"+str(e)
 elif harp==7115:
  r["static_path_fallback_no_cea_metadata"]=True
  r["files"]=[download(harp,s,"/SUM98/D978177825/S00000/"+s+".fits") for s in SEGS]
 RESULT["harps"].append(r);out("harp_"+str(harp)+".json",r)
for date in ("2017-09-03","2017-09-04","2017-09-05","2024-05-13"):
 d=datetime.date.fromisoformat(date)
 for cat,suffix in (("solar_region_summaries","SRS"),("solar_event_reports","events")):
  name=d.strftime("%Y%m%d")+suffix+".txt"
  url=f"https://www.ngdc.noaa.gov/stp/space-weather/swpc-products/daily_reports/{cat}/{d:%Y/%m/}{name}"
  r={"url":url,"state":"SOURCE_UNAVAILABLE"}
  try:
   raw,u=get(url,cap=1000000)
   if len(raw)<50 or b"<html" in raw[:200].lower():raise ValueError("INVALID_NOAA_TEXT")
   folder=O/"noaa";folder.mkdir(exist_ok=True)
   (folder/name).write_bytes(raw)
   r.update(state="SOURCE_BYTES_ACQUIRED",bytes=len(raw),sha256=sha(raw))
  except Exception as e:r["error"]=type(e).__name__+":"+str(e)[:180]
  RESULT["noaa"].append(r)
out("summary.json",RESULT)
print(json.dumps({"harps":[{"harp":h["harp"],"metadata":h["source_metadata_ok"],"fits":sum(f["state"]=="SOURCE_FITS_ACQUIRED" for f in h.get("files",[]))} for h in RESULT["harps"]],"noaa_files":sum(r["state"]=="SOURCE_BYTES_ACQUIRED" for r in RESULT["noaa"]),"science_gate_accepted":False},indent=2))
