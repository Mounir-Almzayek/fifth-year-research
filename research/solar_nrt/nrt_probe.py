#!/usr/bin/env python3
"""Bounded, public-only NRT SHARP CEA source audit. No model acceptance claims."""
import csv,datetime,hashlib,json,math,pathlib,re,sys,time,unittest
import requests

SER="hmi.sharp_cea_720s_nrt"
KEYS="HARPNUM,NOAA_ARS,NOAA_AR,NOAA_NUM,T_REC,T_OBS,QUALITY,CRVAL1,CRVAL2"
HOSTS=("https://jsoc.stanford.edu","https://jsoc1.stanford.edu","http://jsoc.stanford.edu")
SEEDS=[
 (7115,"12673","2017.09.04_00:00:00_TAI",118.228928,-9.243549),
 (7122,"12677","2017.09.04_00:00:00_TAI",49.260674,18.191),
 (9781,"13363","2023.07.16_00:36:00_TAI",339.479919,-21.735651),
 (9802,"13372","2023.07.16_00:36:00_TAI",275.60376,20.586),
 (11170,"13667","2024.05.13_00:36:00_TAI",301.526276,29.6714),
 (11220,"13677","2024.05.13_00:36:00_TAI",313.324921,-3.6445)
]
OUT=pathlib.Path("solar-nrt-evidence")
OUT.mkdir(exist_ok=True)
S=requests.Session();S.headers["User-Agent"]="Solar-Intelligence-scientific-NRT-read-only/1.0"
def save(name,obj):
 (OUT/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,default=str))
def ids(s):
 return {int(v) for v in str(s or "").split(",") if v.strip().isdigit() and int(v)>0}
def distance(a,b,c,d):
 a,b,c,d=map(math.radians,(a,b,c,d))
 v=math.sin(b)*math.sin(d)+math.cos(b)*math.cos(d)*math.cos(a-c)
 return math.degrees(math.acos(max(-1.,min(1.,v))))
def columns(doc,key="keywords"):
 if int(doc.get("status",-1))!=0:raise ValueError("DRMS_NOT_SUCCESS:"+str(doc.get("status")))
 n=int(doc.get("count",0))
 if n<0 or n>400:raise ValueError("UNSAFE_ROW_COUNT")
 rows=[{} for _ in range(n)]
 for col in doc.get(key,[]):
  if len(col.get("values",[]))!=n:raise ValueError("COLUMN_LENGTH_ERROR")
  for i,v in enumerate(col["values"]):rows[i][col["name"]]=v
 return rows
def get(url,params=None,maxbytes=2_500_000):
 with S.get(url,params=params,timeout=(8,20),stream=True) as r:
  r.raise_for_status()
  data=bytearray()
  for block in r.iter_content(131072):
   data.extend(block)
   if len(data)>maxbytes:raise ValueError("SOURCE_SIZE_LIMIT")
  return bytes(data),r.url
def query(epoch):
 dataset=f"{SER}[][{epoch}]"
 params={"op":"rs_list","ds":dataset,"key":KEYS,"seg":"continuum,Br"}
 attempts=[]
 for host in HOSTS:
  try:
   raw,url=get(host+"/cgi-bin/ajax/jsoc_info",params)
   doc=json.loads(raw)
   kw=columns(doc)
   segments=columns(doc,"segments")
   if len(kw)!=len(segments):raise ValueError("SEGMENT_COUNT_MISMATCH")
   return {"ok":True,"url":url,"sha256":hashlib.sha256(raw).hexdigest(),
           "data":doc,"rows":kw,"segments":segments,"attempts":attempts}
  except Exception as e:
   attempts.append({"host":host,"error":type(e).__name__+":"+str(e)[:230]})
 return {"ok":False,"attempts":attempts,"dataset":dataset}
def match(seed,rows):
 definitive_harp,noaa,epoch,lon,lat=seed
 possible=[]
 for index,row in enumerate(rows):
  try:
   gap=distance(lon,lat,float(row["CRVAL1"]),float(row["CRVAL2"]))
   rel=bool(ids(noaa)&ids(row.get("NOAA_ARS")))
   possible.append({"index":index,"nrt_harpnum":row["HARPNUM"],"nrt_noaa":row.get("NOAA_ARS",""),
                    "angular_gap_deg":round(gap,4),"noaa_overlap":rel,
                    "quality":row.get("QUALITY"),"source_t_rec":row.get("T_REC"),
                    "source_t_obs":row.get("T_OBS")})
  except Exception as e:continue
 possible.sort(key=lambda v:(not v["noaa_overlap"],v["angular_gap_deg"]))
 close=[x for x in possible if x["angular_gap_deg"]<=10.]
 noaa_near=[x for x in close if x["noaa_overlap"]]
 state=("NO_POSITION_MATCH" if not close else "UNIQUE_NOAA_AND_POSITION_CANDIDATE" if len(noaa_near)==1
        else "UNIQUE_GEOMETRIC_CANDIDATE" if len(close)==1 else "AMBIGUOUS_SOURCE_MATCH")
 # Candidates are never automatically physical-family accepted.
 return {"discovery_status":state,"plausible_count":len(close),"candidates":possible[:30],
         "eligible_as_issued":False,"historical_publication_time_utc":"UNKNOWN"}
def fetch_segment(path,filename):
 if not isinstance(path,str) or not re.fullmatch(r"/SUM[0-9]+/D[0-9]+/S[0-9]+/[A-Za-z0-9_]+\.fits",path):
  return {"state":"UNSAFE_OR_MISSING_SOURCE_PATH"}
 attempts=[]
 for host in ("https://jsoc1.stanford.edu","https://jsoc.stanford.edu"):
  try:
   data,url=get(host+path,maxbytes=12_000_000)
   if len(data)%2880 or not data.startswith(b"SIMPLE  ="):
    raise ValueError("BAD_FITS_STRUCTURE")
   filename.parent.mkdir(parents=True,exist_ok=True)
   filename.write_bytes(data)
   info={"state":"NRT_SOURCE_FITS_ACQUIRED","sha256":hashlib.sha256(data).hexdigest(),
         "bytes":len(data),"source_url":url}
   try:
    import astropy.io.fits as fits,numpy as np,io
    with fits.open(io.BytesIO(data),checksum=True,memmap=False) as hdus:
     imgs=[np.asarray(h.data) for h in hdus if getattr(h.data,"ndim",0)==2]
     info["shape"]=[list(a.shape) for a in imgs]
     info["checksum_verified"]=[h.verify_checksum() if "CHECKSUM" in h.header else None for h in hdus]
   except Exception as e:info["pixel_issue"]=type(e).__name__+":"+str(e)[:160]
   return info
  except Exception as e:attempts.append(type(e).__name__+":"+str(e)[:190])
 return {"state":"NRT_SEGMENT_FETCH_FAILED","attempts":attempts}
class LogicTests(unittest.TestCase):
 def test_wrap_longitude(self):self.assertAlmostEqual(distance(359,0,1,0),2,places=4)
 def test_exact_position(self):self.assertAlmostEqual(distance(10,10,10,10),0,places=4)
 def test_noaa_ids(self):self.assertEqual(ids("12673,12674"),{12673,12674})
 def test_missing_ids(self):self.assertEqual(ids("MISSING"),set())
 def test_no_source_rows(self):self.assertEqual(match(SEEDS[0],[])["discovery_status"],"NO_POSITION_MATCH")
 def test_no_definitive_id_inheritance(self):
  row={"HARPNUM":"9999","NOAA_ARS":"12673","CRVAL1":"118.2","CRVAL2":"-9.24"}
  self.assertEqual(match(SEEDS[0],[row])["candidates"][0]["nrt_harpnum"],"9999")
 def test_no_operational_release_claim(self):
  row={"HARPNUM":"9999","NOAA_ARS":"12673","CRVAL1":"118.2","CRVAL2":"-9.24"}
  self.assertFalse(match(SEEDS[0],[row])["eligible_as_issued"])
 def test_bad_col_length(self):
  with self.assertRaises(ValueError):columns({"status":0,"count":2,"keywords":[{"name":"x","values":["1"]}]})
 def test_bad_drms_status(self):
  with self.assertRaises(ValueError):columns({"status":7,"count":0})
 def test_good_drms_columns(self):
  self.assertEqual(columns({"status":0,"count":1,"keywords":[{"name":"HARPNUM","values":["5"]}]}),[{"HARPNUM":"5"}])
if __name__=="__main__":
 suite=unittest.defaultTestLoader.loadTestsFromTestCase(LogicTests)
 result=unittest.TextTestRunner(verbosity=2).run(suite)
 if not result.wasSuccessful():raise SystemExit(1)
 status={"science_gate":"G1_IN_PROGRESS_G2_G3_G4_NOT_ACCEPTED",
         "as_issued_proven":False,"source_series":SER,
         "query_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
         "source_epochs":[],"sample_matches":[]}
 by_epoch={}
 for epoch in dict.fromkeys(x[2] for x in SEEDS):
  response=query(epoch)
  status["source_epochs"].append({"t_rec_tai":epoch,"queried":response["ok"],
                                "url":response.get("url"),
                                "attempts":response["attempts"],
                                "source_sha256":response.get("sha256"),
                                "record_count":len(response["rows"]) if response["ok"] else None})
  if response["ok"]:
   save("jsoc_source_"+epoch.replace(".","_").replace(":","")+".json",response["data"])
  by_epoch[epoch]=response
 downloaded=0
 for seed in SEEDS:
  resp=by_epoch[seed[2]]
  rec={"definitive_harpnum":seed[0],"definitive_noaa":seed[1],
       "t_rec":seed[2],"source_class":"DEFINITIVE_ONLY",
       "nrt_publication_utc":"UNKNOWN",
       "eligible_as_issued":False}
  if resp["ok"]:
   m=match(seed,resp["rows"])
   rec.update(m)
   c=[x for x in m["candidates"] if x["angular_gap_deg"]<=10 and x["noaa_overlap"]]
   if len(c)==1 and downloaded<2:
    nr=resp["segments"][c[0]["index"]]
    rec["source_segments"]={}
    for seg in ("continuum","Br"):
     fname=OUT/f"nrt_{c[0]['nrt_harpnum']}_{seed[0]}"/f"{seg}.fits"
     rec["source_segments"][seg]=fetch_segment(nr.get(seg),fname)
    downloaded+=1
  else:
   rec["discovery_status"]="PUBLIC_JSOC_QUERY_UNAVAILABLE"
  status["sample_matches"].append(rec)
 save("nrt_g1_source_audit.json",status)
 with (OUT/"nrt_g1_candidates.csv").open("w",newline="",encoding="utf-8") as f:
  w=csv.writer(f);w.writerow(["def_harp","def_noaa","t_rec","state","nrt_harp","nrt_noaa","great_circle_deg","noaa_overlap","as_issued_verified"])
  for rec in status["sample_matches"]:
   for c in (rec.get("candidates") or [{}]):
    w.writerow([rec["definitive_harpnum"],rec["definitive_noaa"],rec["t_rec"],rec["discovery_status"],
               c.get("nrt_harpnum",""),c.get("nrt_noaa",""),c.get("angular_gap_deg",""),c.get("noaa_overlap",""),False])
 print(json.dumps({"science_gate":status["science_gate"],"source_epochs":status["source_epochs"],
                   "candidate_statuses":[[r["definitive_harpnum"],r["discovery_status"]] for r in status["sample_matches"]]},indent=2))


# G1_WCS_DATE_V2_APPEND — 2026-11-21: original-source WCS and DATE audit.
# Bounded to two pairs, all original bytes/provenance preserved; never G1/G3 acceptance.
def wcs_date_v2():
 import io,os,numpy as np
 from astropy.io import fits
 from astropy.wcs import WCS
 out=pathlib.Path("solar-nrt-evidence")/"wcs_date_v2"
 out.mkdir(parents=True,exist_ok=True)
 pairs=[("2017.09.04_00:00:00_TAI",5634,7115),
        ("2023.07.16_00:36:00_TAI",8606,9781)]
 groups=["HARPNUM,NOAA_AR,NOAA_NUM,NOAA_ARS,T_REC,T_OBS,QUALITY,CRPIX1,CRPIX2,CRVAL1,CRVAL2,CDELT1,CDELT2,CTYPE1,CTYPE2",
         "DATE,DATE_S,DATE_B,DATE__OBS,DATE_OBS",
         "RSUN_REF,RSUN_OBS,DSUN_OBS,CUNIT1,CUNIT2,WCSNAME,PV2_1"]
 def qry(series,harp,epoch,keys,segs=None):
  ds=f"{series}[{harp}][{epoch}]"
  params={"op":"rs_list","ds":ds,"key":keys}
  if segs:params["seg"]=segs
  attempts=[]
  for host in HOSTS:
   try:
    raw,url=get(host+"/cgi-bin/ajax/jsoc_info",params,maxbytes=4_000_000)
    obj=json.loads(raw)
    if obj.get("status")!=0:raise ValueError("DRMS_STATUS_"+str(obj.get("status"))+":"+str(obj.get("error",""))[:120])
    row=columns(obj)
    if len(row)!=1:raise ValueError("EXPECTED_ONE_ROW_GOT_"+str(len(row)))
    segrows=columns(obj,"segments") if segs else []
    return {"ok":True,"source_url":url,"source_sha256":hashlib.sha256(raw).hexdigest(),
            "original_drms_json":obj,"row":row[0],
            "segments":segrows[0] if len(segrows)==1 else {},"attempts":attempts}
   except Exception as ex:attempts.append({"host":host,"error":type(ex).__name__+":"+str(ex)[:180]})
  return {"ok":False,"attempts":attempts,"dataset":ds,"requested_keys":keys}
 def source_wcs(d,shape):
  # All metadata comes from raw DRMS query for the exact source record; never
  # substitute reference values from the other series.
  keys=["CRPIX1","CRPIX2","CRVAL1","CRVAL2","CDELT1","CDELT2","CTYPE1","CTYPE2"]
  if not all(k in d for k in keys):raise ValueError("WCS_KEYS_NOT_AVAILABLE")
  h=fits.Header()
  h["NAXIS"]=2;h["NAXIS1"]=int(shape[1]);h["NAXIS2"]=int(shape[0])
  for k in keys:
   val=d[k]
   h[k]=float(val) if k.startswith(("CRPIX","CRVAL","CDELT")) else str(val)
  h["CUNIT1"]="deg";h["CUNIT2"]="deg"
  return WCS(h,fix=False,relax=True)
 def cutout_data(b):
  with fits.open(io.BytesIO(b),checksum=True,memmap=False) as hdus:
   arrays=[np.asarray(h.data,dtype=np.float64) for h in hdus if getattr(h.data,"ndim",0)==2]
   if len(arrays)!=1:raise ValueError("EXPECTED_ONE_IMAGE_HDU")
   valid=all(h.verify_checksum()==1 for h in hdus if "CHECKSUM" in h.header)
   return arrays[0],valid
 def bilinear(img,x,y):
  a=np.floor(x).astype(int); b=np.floor(y).astype(int)
  xdiff=x-a;ydiff=y-b
  return ((1-xdiff)*(1-ydiff)*img[b,a]+xdiff*(1-ydiff)*img[b,a+1]
          +(1-xdiff)*ydiff*img[b+1,a]+xdiff*ydiff*img[b+1,a+1])
 results=[]
 for epoch,nrt_harp,def_harp in pairs:
  rec={"t_rec_tai":epoch,"nrt_harp":nrt_harp,"definitive_harp":def_harp,
       "publication_utc_verified":False,"science_gate":"G1_IN_PROGRESS"}
  meta={}
  for tag,series,h in (("nrt","hmi.sharp_cea_720s_nrt",nrt_harp),
                       ("definitive","hmi.sharp_cea_720s",def_harp)):
   a=qry(series,h,epoch,groups[0],"continuum,Br")
   b=qry(series,h,epoch,groups[1])
   c=qry(series,h,epoch,groups[2])
   meta[tag]={"base":a,"date_fields":b,"auxiliary":c}
   save("wcs_date_v2/source_"+tag+"_"+str(h)+".json",meta[tag])
  rec["query_status"]={k:{kind:v[k]["ok"] for kind in v} for k,v in meta.items()}
  rec["core_fields"]={k:v["base"].get("row",{}) for k,v in meta.items()}
  rec["date_fields"]={k:v["date_fields"].get("row",{}) for k,v in meta.items()}
  rec["DATE_processing_not_publication"]=True
  try:
   if not all(v["base"]["ok"] for v in meta.values()):raise ValueError("CORE_SOURCE_QUERY_FAILED")
   for segment in ("continuum","Br"):
    objects={}
    for tag,v in meta.items():
     relative=v["base"]["segments"].get(segment)
     if not isinstance(relative,str) or not re.fullmatch(r"/SUM[0-9]+/D[0-9]+/S[0-9]+/[A-Za-z0-9_]+\.fits",relative):
      raise ValueError("UNSAFE_OR_MISSING_SEGMENT_"+tag+"_"+segment)
     errors=[]
     for host in HOSTS:
      try:
       binary,url=get(host+relative,maxbytes=9_000_000)
       if len(binary)%2880!=0 or not binary.startswith(b"SIMPLE  ="):raise ValueError("NOT_FITS")
       arr,chksum=cutout_data(binary)
       objects[tag]={"array":arr,"checksum":chksum,"source_url":url,
         "sha256":hashlib.sha256(binary).hexdigest(),"shape":list(arr.shape)}
       break
      except Exception as ex:errors.append(str(ex)[:120])
     if tag not in objects:raise ValueError("SOURCE_DOWNLOAD_FAILURE_"+tag+"_"+str(errors))
    nd=objects["nrt"];df=objects["definitive"]
    w1=source_wcs(meta["definitive"]["base"]["row"],df["array"].shape)
    w2=source_wcs(meta["nrt"]["base"]["row"],nd["array"].shape)
    ny,nx=df["array"].shape
    step=max(1,int(np.sqrt(nx*ny/25000.)))
    ys,xs=np.mgrid[0:ny:step,0:nx:step]
    xx=xs.astype(float).ravel();yy=ys.astype(float).ravel()
    lon,lat=w1.all_pix2world(xx,yy,0)
    x2,y2=w2.all_world2pix(lon,lat,0)
    lon2,lat2=w2.all_pix2world(x2,y2,0)
    wrap=lambda a:((a+180)%360)-180
    deviation=np.hypot(wrap(lon2-lon)*np.cos(np.deg2rad(lat)),lat2-lat)
    valid=np.isfinite(x2)&np.isfinite(y2)&np.isfinite(deviation)&(deviation<1.e-4)
    valid&=(x2>=0)&(y2>=0)&(x2<nd["array"].shape[1]-1)&(y2<nd["array"].shape[0]-1)
    overlap=int(valid.sum())
    metrics={"source_files":{k:{"source_url":o["source_url"],"sha256":o["sha256"],
           "shape":o["shape"],"checksum_valid":o["checksum"]} for k,o in objects.items()},
        "attempted_grid_samples":len(xx),"mapped_to_both_pixels":overlap,
        "overlap_fraction":overlap/max(1,len(xx)),
        "wcs_roundtrip_max_degrees":float(np.nanmax(deviation[valid])) if overlap else None,
        "method":"Exact-record DRMS CEA WCS; astropy/WCSLIB pix_to_world/world_to_pix; bilinear on matched coordinates; nonindependent sample pixels"}
    if overlap>=100:
     sel=np.where(valid)[0]
     in_def=df["array"][yy[sel].astype(int),xx[sel].astype(int)]
     in_nrt=bilinear(nd["array"],x2[sel],y2[sel])
     good=np.isfinite(in_def)&np.isfinite(in_nrt)
     if good.sum()>=100:
      one=in_def[good];two=in_nrt[good]
      metrics.update({"valid_co_registered_pixels":int(good.sum()),"pearson":float(np.corrcoef(one,two)[0,1]),
        "median_difference_def_minus_nrt":float(np.median(one-two)),
        "median_abs_difference":float(np.median(np.abs(one-two))),
        "definite_median":float(np.median(one)),"nrt_median":float(np.median(two))})
    rec.setdefault("co_registration",{})[segment]=metrics
  except Exception as ex:
   rec["co_registration_blocker"]=type(ex).__name__+":"+str(ex)[:700]
  results.append(rec)
 save("wcs_date_v2/summary.json",{"source":"ORIGINAL_JSOC_EXACT_RECORDS","results":results,
   "nrt_publication_utc_verified":False,"scientific_gates":{"G1":"IN_PROGRESS","G2":"NOT_ACCEPTED","G3":"NOT_ACCEPTED","G4":"NOT_ACCEPTED"}})
 print("G1 WCS DATE v2 summary:",json.dumps([{"pair":f"{r['definitive_harp']}/{r['nrt_harp']}",
         "meta":r["query_status"],"registered":list(r.get("co_registration",{})),
         "blocker":r.get("co_registration_blocker")} for r in results]))
try:
 wcs_date_v2()
except Exception as _g1_v2_exception:
 save("wcs_date_v2/fatal_error.json",{"error":type(_g1_v2_exception).__name__+":"+str(_g1_v2_exception)})
 print("G1 WCS v2 failed safely:",_g1_v2_exception)
