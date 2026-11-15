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
 for host in ("http://jsoc.stanford.edu","https://jsoc.stanford.edu","https://jsoc1.stanford.edu"):
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


# Export scientific pixel arrays of the *actual* original 2017 FITS in the
# exact physical values decoded by Astropy, to permit CV QA in offline environments.
try:
 import numpy as np
 from astropy.io import fits
 dir7115=O/"fits_7115"
 sources={s:dir7115/(s+".fits") for s in SEGS}
 if all(path.exists() for path in sources.values()):
  raw_arrays={}
  for seg,path in sources.items():
   with fits.open(path,memmap=False) as hdus:
    a=[np.asarray(x.data) for x in hdus if getattr(x.data,"ndim",0)==2]
    if len(a)!=1:raise ValueError("EXPECTED_ONE_2D_HDU_"+seg)
    raw_arrays[seg]=a[0].astype("float32")
  np.savez_compressed(O/"hmi_7115_original_pixel_arrays.npz",**raw_arrays)
  out("decoded_array_manifest.json",{"origin":"JSOC1 original source HARP7115 2017-09-04", "note":"NUMPY ARRAYS NOT INDEPENDENT GROUND TRUTH","arrays":{k:{"shape":list(v.shape),"min":float(np.nanmin(v)),"max":float(np.nanmax(v)),"finite_fraction":float(np.isfinite(v).mean())} for k,v in raw_arrays.items()}})
except Exception as e:
 out("decoded_array_error.json",{"error":type(e).__name__+": "+str(e)})

# Auxiliary 2014 optical test case, obtained from third-party GitHub; the
# published mask is algorithm-generated and explicitly NOT independent truth.
try:
 from io import BytesIO
 import numpy as np
 from astropy.io import fits
 root="https://raw.githubusercontent.com/wangzai15/sunspot-mat/main/"
 paths={"continuum":"norm_sub_hmi.ic_45s.20140811_083000_TAI.2.continuum.fits","auto_mask":"ternary_mask.npy"}
 sample={"role":"THIRD_PARTY_CV_ONLY","independent_reference_mask":False,"files":{}}
 for kind,path in paths.items():
  b,u=get(root+path,cap=7000000)
  info={"source_url":u,"bytes":len(b),"sha256":sha(b)}
  if kind=="continuum":
   with fits.open(BytesIO(b),memmap=False) as hdus:
    img=np.asarray(hdus[0].data)
    info.update(shape=list(img.shape),finite_fraction=float(np.isfinite(img).mean()),
                min=float(np.nanmin(img)),max=float(np.nanmax(img)))
  else:
   mask=np.load(BytesIO(b),allow_pickle=False)
   unique,counts=np.unique(mask,return_counts=True)
   info.update(shape=list(mask.shape),class_histogram={str(k):int(v) for k,v in zip(unique,counts)})
  sample["files"][kind]=info
 out("thirdparty_2014_source_audit.json",sample)
except Exception as e:
 out("thirdparty_2014_source_error.json",{"type":type(e).__name__,"error":str(e)[:400]})


# Separate 2024 multi-NOAA HARP stress case, independently documented on
# the project's 2026-11-06 original JSOC acquisition. Not an eligible
# one-NOAA flare-label row; valuable for real optical/magnetic CV QA.
T24="2024.05.10_00:00:00_TAI";H24=11149
q24=query(f"{SER}[{H24}][{T24}]",FIELDS,",".join(SEGS))
out("query_11149_2024.json",q24)
test24={"harp":H24,"t_rec":T24,"status":"MULTI_NOAA_CV_STRESS_ONLY","metadata_retrieved":q24["ok"]}
if q24["ok"]:
 try:
  rr=rows(q24["data"],"keywords");ss=rows(q24["data"],"segments")
  test24["records"]=rr
  if len(rr)==len(ss)==1:test24["files"]=[download(H24,s,ss[0].get(s)) for s in SEGS]
  else:test24["problem"]="NOT_SINGLE_CEA_RECORD"
 except Exception as e:test24["problem"]=type(e).__name__+": "+str(e)
else:
 test24["static_path_from_earlier_verified_source"]=True
 test24["files"]=[download(H24,s,"/SUM97/D1753927374/S00000/"+s+".fits") for s in SEGS]
out("harp_11149_2024.json",test24)
try:
 from astropy.io import fits
 import numpy as np
 dp=O/"fits_11149"
 if all((dp/(s+".fits")).exists() for s in SEGS):
  data24={}
  for seg in SEGS:
   with fits.open(dp/(seg+".fits"),memmap=False) as hdus:
    arrays=[np.asarray(h.data) for h in hdus if getattr(h.data,"ndim",0)==2]
    if len(arrays)!=1:raise ValueError("ONE_IMAGE_HDU_REQUIRED_"+seg)
    data24[seg]=arrays[0].astype("float32")
  np.savez_compressed(O/"hmi_11149_2024_original_pixel_arrays.npz",**data24)
  out("decoded_11149_2024_manifest.json",{"scientific_label":"UNKNOWN_MULTI_NOAA","arrays":{s:{"shape":list(a.shape),"finite_fraction":float(np.isfinite(a).mean())} for s,a in data24.items()}})
except Exception as e:out("decoded_11149_2024_error.json",{"error":str(e)})

# Additional original daily NOAA source bytes to cross-check the May 2024
# event study, without implying complete negative-label coverage.
for date in ("2024-05-08","2024-05-09","2024-05-10","2024-05-11","2024-05-12"):
 d=datetime.date.fromisoformat(date)
 for cat,suffix in (("solar_region_summaries","SRS"),("solar_event_reports","events")):
  name=d.strftime("%Y%m%d")+suffix+".txt"
  url=f"https://www.ngdc.noaa.gov/stp/space-weather/swpc-products/daily_reports/{cat}/{d:%Y/%m/}{name}"
  entry={"url":url,"status":"SOURCE_UNAVAILABLE"}
  try:
   data,u=get(url,cap=1000000)
   if len(data)<50 or b"<html" in data[:200].lower():raise ValueError("INVALID_SOURCE_TEXT")
   folder=O/"noaa";folder.mkdir(exist_ok=True);(folder/name).write_bytes(data)
   entry.update(status="SOURCE_BYTES_ACQUIRED",bytes=len(data),sha256=sha(data))
  except Exception as e:entry["error"]=type(e).__name__+": "+str(e)[:180]
  out("noaa_may_2024_"+name+".json",entry)


# Cross-epoch, deliberately bounded official CEA candidate discovery.
# Candidate NOAA IDs are sourced from previously archived NOAA SRS dates,
# NOT asserted to be independently verified until the full original join.
CROSS_EPOCH=[
 ("2024.05.13_00:36:00_TAI",13674,"may_2024_target_M1_candidate"),
 ("2023.07.16_00:36:00_TAI",13371,"jul_2023_control_candidate"),
 ("2025.03.22_00:36:00_TAI",14030,"mar_2025_control_candidate"),
]
multi_epoch=[]
for epoch,noaa,title in CROSS_EPOCH:
 rec={"epoch":epoch,"noaa_candidate":noaa,"title":title,"scientifically_accepted":False,"files":[]}
 rev=query(f"{SER}[][{epoch}]","HARPNUM,NOAA_AR,NOAA_NUM,NOAA_ARS,T_REC,QUALITY")
 out(f"epoch_reverse_{title}.json",rev)
 if not rev["ok"]:
  rec["state"]="SOURCE_UNAVAILABLE";rec["errors"]=rev["attempts"];multi_epoch.append(rec);continue
 try:
  raw_rows=rows(rev["data"],"keywords")
  def any_noaa(row):
   return noaa in [int(x) for x in str(row.get("NOAA_ARS","")).split(",") if x.strip().isdigit()]
  matches=[row for row in raw_rows if any_noaa(row)]
  rec["matching_reverse_rows"]=matches
  if len(matches)!=1 or int(matches[0].get("NOAA_NUM",0))!=1 or int(matches[0].get("NOAA_AR",0))!=noaa:
   rec["state"]="NO_UNIQUE_ONE_NOAA_HARP";multi_epoch.append(rec);continue
  harp=int(matches[0]["HARPNUM"])
  rec["harpnum"]=harp
  exact=query(f"{SER}[{harp}][{epoch}]",FIELDS,",".join(SEGS))
  out(f"epoch_exact_{title}.json",exact)
  if not exact["ok"]:
   rec["state"]="EXACT_SOURCE_UNAVAILABLE";multi_epoch.append(rec);continue
  er=rows(exact["data"],"keywords");sg=rows(exact["data"],"segments")
  if len(er)!=1 or len(sg)!=1 or int(er[0]["NOAA_NUM"])!=1 or int(er[0]["NOAA_AR"])!=noaa:
   rec["state"]="EXACT_JOIN_REJECTED";multi_epoch.append(rec);continue
  rec["record"]=er[0]
  rec["files"]=[download(harp,s,sg[0].get(s)) for s in SEGS]
  rec["state"]="SOURCE_ACQUIRED_NEEDS_SCIENCE_REVIEW" if all(x["state"]=="SOURCE_FITS_ACQUIRED" for x in rec["files"]) else "INCOMPLETE_SEGMENTS"
  if rec["state"]=="SOURCE_ACQUIRED_NEEDS_SCIENCE_REVIEW":
   try:
    import numpy as np
    from astropy.io import fits
    src=O/("fits_"+str(harp))
    dec={}
    for seg in SEGS:
     with fits.open(src/(seg+".fits"),memmap=False) as hdus:
      candidate=[np.asarray(h.data) for h in hdus if getattr(h.data,"ndim",0)==2]
      if len(candidate)!=1:raise ValueError("ONE_IMAGE_REQUIRED_"+seg)
      dec[seg]=candidate[0].astype("float32")
    np.savez_compressed(O/("hmi_"+str(harp)+"_epoch_pixel_arrays.npz"),**dec)
    rec["decoded"]={k:{"shape":list(v.shape),"finite_fraction":float(np.isfinite(v).mean())} for k,v in dec.items()}
   except Exception as e:rec["decode_error"]=type(e).__name__+":"+str(e)[:220]
 except Exception as e:
  rec["state"]="AUDIT_ERROR";rec["error"]=type(e).__name__+":"+str(e)[:220]
 multi_epoch.append(rec)
out("cross_epoch_discovery.json",multi_epoch)

# Three-day source-horizon review files only; no unearned negative-label
# acceptances. NOAA SRS and edited XRA must be adjudicated as-of separately.
for d in ("2024-05-14","2024-05-15","2023-07-15","2023-07-16","2023-07-17","2025-03-21","2025-03-22","2025-03-23"):
 day=datetime.date.fromisoformat(d)
 for category,sfx in (("solar_region_summaries","SRS"),("solar_event_reports","events")):
  name=day.strftime("%Y%m%d")+sfx+".txt"
  url=f"https://www.ngdc.noaa.gov/stp/space-weather/swpc-products/daily_reports/{category}/{day:%Y/%m/}{name}"
  state={"url":url,"status":"SOURCE_UNAVAILABLE"}
  try:
   contents,got=get(url,cap=1000000)
   if len(contents)<50 or b"<html" in contents[:200].lower():raise ValueError("SOURCE_TEXT_INVALID")
   (O/"noaa").mkdir(exist_ok=True);(O/"noaa"/name).write_bytes(contents)
   state.update(status="SOURCE_BYTES_ACQUIRED",bytes=len(contents),sha256=sha(contents))
  except Exception as e:state["error"]=type(e).__name__+":"+str(e)[:180]
  out("noaa_extended_"+name+".json",state)


# Deterministically selected cross-epoch one-NOAA CEA source candidates, from
# actual 2023/2024/2025 same-time official CEA reverse query. Mixed complexity
# and epoch; not selected by future flare labels, no approved family graph yet.
CLEAN_SAMPLE=[
 ("2023.07.16_00:36:00_TAI",9781,13363),
 ("2023.07.16_00:36:00_TAI",9802,13372),
 ("2024.05.13_00:36:00_TAI",11170,13667),
 ("2024.05.13_00:36:00_TAI",11220,13677),
 ("2025.03.22_00:36:00_TAI",12906,14026),
 ("2025.03.22_00:36:00_TAI",12941,14035),
]
clean_samples=[]
for epoch,harp,noaa in CLEAN_SAMPLE:
 result={"epoch":epoch,"harp":harp,"noaa":noaa,"state":"UNKNOWN","accepted":False}
 reverse=query(f"{SER}[][{epoch}]","HARPNUM,NOAA_AR,NOAA_ARS,NOAA_NUM,T_REC,QUALITY")
 if not reverse["ok"]:
  result["state"]="SOURCE_UNAVAILABLE";clean_samples.append(result);continue
 rr=rows(reverse["data"],"keywords")
 valid=[r for r in rr if r.get("NOAA_AR")==str(noaa) and r.get("NOAA_ARS")==str(noaa) and r.get("NOAA_NUM")=="1" and r.get("HARPNUM")==str(harp)]
 all_matches=[r for r in rr if str(noaa) in str(r.get("NOAA_ARS","")).split(",")]
 if len(valid)!=1 or len(all_matches)!=1:
  result["state"]="REVERSE_IDENTITY_AMBIGUOUS";clean_samples.append(result);continue
 exact=query(f"{SER}[{harp}][{epoch}]",FIELDS,",".join(SEGS))
 out(f"clean_exact_{harp}.json",exact)
 if not exact["ok"]:
  result["state"]="EXACT_CEA_SOURCE_UNAVAILABLE";clean_samples.append(result);continue
 try:
  records=rows(exact["data"],"keywords");segments=rows(exact["data"],"segments")
  if len(records)!=1 or len(segments)!=1 or records[0].get("NOAA_NUM")!="1" or records[0].get("NOAA_ARS")!=str(noaa):
   raise ValueError("EXACT_CEA_IDENTITY_REJECTED")
  result["record"]=records[0]
  result["files"]=[download(harp,s,segments[0].get(s)) for s in SEGS]
  if not all(x["state"]=="SOURCE_FITS_ACQUIRED" for x in result["files"]):
   result["state"]="INCOMPLETE_SEGMENTS"
  else:
   result["state"]="SOURCE_ACQUIRED_NEEDS_PHYSICAL_REVIEW"
   import numpy as np
   from astropy.io import fits
   dec={}
   for seg in SEGS:
    with fits.open(O/f"fits_{harp}"/f"{seg}.fits",memmap=False) as hdus:
     candidates=[np.asarray(h.data) for h in hdus if getattr(h.data,"ndim",0)==2]
     if len(candidates)!=1:raise ValueError("SOURCE_IMAGE_HDU_REQUIRED")
     dec[seg]=candidates[0].astype("float32")
   np.savez_compressed(O/f"hmi_{harp}_clean_candidate_pixels.npz",**dec)
   result["decoded"]={k:{"shape":list(v.shape),"finite_fraction":float(np.isfinite(v).mean())} for k,v in dec.items()}
 except Exception as e:
  result["state"]="SOURCE_OR_DECODE_QA_ERROR"
  result["error"]=type(e).__name__+":"+str(e)[:260]
 clean_samples.append(result)
out("clean_cross_epoch_sources.json",clean_samples)


# Close remaining original HARP7122 raw-array decoding evidence and obtain
# the 2024 multi-NOAA reverse metadata independently (never accepted labels).
try:
 import numpy as np
 from astropy.io import fits
 target=O/"fits_7122"
 pixel={}
 if all((target/(s+".fits")).is_file() for s in SEGS):
  for seg in SEGS:
   with fits.open(target/(seg+".fits"),memmap=False) as hdus:
    dd=[np.asarray(h.data) for h in hdus if getattr(h.data,"ndim",0)==2]
    if len(dd)!=1:raise ValueError("HARP7122_EXPECTED_ONE_2D_HDU")
    pixel[seg]=dd[0].astype("float32")
  np.savez_compressed(O/"hmi_7122_original_pixel_arrays.npz",**pixel)
  out("decoded_7122_manifest.json",{"harp":7122,"noaa":12677,"independent_ground_truth":False,
    "segments":{k:{"shape":list(v.shape),"finite_fraction":float(np.isfinite(v).mean())} for k,v in pixel.items()}})
except Exception as e:
 out("decoded_7122_error.json",{"error":type(e).__name__+":"+str(e)[:220]})
reverse2024=query(f"{SER}[][2024.05.10_00:00:00_TAI]","HARPNUM,NOAA_AR,NOAA_ARS,NOAA_NUM,T_REC,QUALITY")
out("reverse_cea_2024_11149.json",reverse2024)


# Official SHARP scalar-keyword descriptors from the exact same definitive
# CEA records. Source values are instrument's precomputed magnetic baseline;
# do not backpropagate label knowledge or claim original optical mask truth.
keyword_keys="HARPNUM,NOAA_AR,NOAA_NUM,NOAA_ARS,T_REC,QUALITY,USFLUX,MEANGAM,MEANGBT,MEANGBH,MEANGBZ,MEANJZD,TOTUSJZ,MEANALP,MEANJZH,TOTUSJH,SAVNCPP,MEANPOT,TOTPOT,MEANSHR,SHRGT45,R_VALUE,AREA_ACR"
archive_records=[
 (7115,"2017.09.04_00:00:00_TAI"),
 (7122,"2017.09.04_00:00:00_TAI"),
 (9781,"2023.07.16_00:36:00_TAI"),
 (9802,"2023.07.16_00:36:00_TAI"),
 (11149,"2024.05.10_00:00:00_TAI"),
 (11170,"2024.05.13_00:36:00_TAI"),
 (11220,"2024.05.13_00:36:00_TAI"),
 (12941,"2025.03.22_00:36:00_TAI"),
]
derived=[]
for harp,epoch in archive_records:
 rec={"harpnum":harp,"t_rec":epoch,"status":"SOURCE_UNAVAILABLE","scientifically_accepted":False}
 r=query(f"{SER}[{harp}][{epoch}]",keyword_keys)
 out(f"sharp_derived_keyword_record_{harp}.json",r)
 if r["ok"]:
  try:
   v=rows(r["data"],"keywords")
   if len(v)==1 and v[0]["HARPNUM"]==str(harp):
    rec.update(status="ACTUAL_JSOC_SOURCE_KEYWORDS",source_sha256=r["sha256"],keywords=v[0])
   else:rec["status"]="UNEXPECTED_RECORD_COUNT"
  except Exception as e:rec["error"]=type(e).__name__+":"+str(e)[:150]
 else:rec["failures"]=r["attempts"]
 derived.append(rec)
out("original_sharp_derived_features.json",derived)
