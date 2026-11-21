#!/usr/bin/env python3
"""Source-byte GOES-16/18 1min XRS source recovery for HARP12941 window.
Never marks negative without instrument quality/continuous coverage review.
"""
import datetime,hashlib,json,pathlib,unittest,re
import requests
OUT=pathlib.Path("solar-goes-xrs-20250322-evidence");OUT.mkdir(exist_ok=True)
BASE="https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes"
S=requests.Session();S.headers["User-Agent"]="Solar-Intelligence-XRS-Source-Custody/1.0"
TARGETS=[(sat,day,f"sci_xrsf-l2-avg1m_g{sat}_d202503{day}_v2-2-0.nc")
 for sat in (16,18) for day in ("22","23")]
def get(url,max_bytes=6_000_000):
 with S.get(url,timeout=(12,90),stream=True) as r:
  r.raise_for_status()
  d=bytearray()
  for piece in r.iter_content(131072):
   d.extend(piece)
   if len(d)>max_bytes:raise ValueError("MAX_BYTES_EXCEEDED")
  return bytes(d)
def inspect(p):
 import netCDF4
 with netCDF4.Dataset(str(p)) as ds:
  variables={}
  for name,var in ds.variables.items():
   variables[name]={"dtype":str(var.dtype),"dimensions":list(var.dimensions),
                    "shape":list(var.shape),"units":str(getattr(var,"units","")),
                    "long_name":str(getattr(var,"long_name",""))[:200],
                    "_FillValue":str(getattr(var,"_FillValue",None)),
                    "flag_values":str(getattr(var,"flag_values",""))[:200],
                    "flag_meanings":str(getattr(var,"flag_meanings",""))[:300]}
  return {"global_attrs":{k:str(ds.getncattr(k))[:450] for k in ds.ncattrs()
                           if k in ("title","source","summary","time_coverage_start","time_coverage_end","processing_level","history","algorithm_version")},
          "dimensions":{k:len(v) for k,v in ds.dimensions.items()},"variables":variables}
class Guard(unittest.TestCase):
 def test_targets(self):self.assertEqual(len(TARGETS),4)
 def test_satellites(self):self.assertEqual({s for s,_,_ in TARGETS},{16,18})
 def test_dates(self):self.assertEqual({d for _,d,_ in TARGETS},{"22","23"})
 def test_filename(self):self.assertTrue(TARGETS[0][2].endswith(".nc"))
if __name__=="__main__":
 t=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Guard))
 if not t.wasSuccessful():raise SystemExit(1)
 out={"source":"NOAA NCEI GOES XRS L2 science-quality 1minute average; original NetCDF",
      "window_utc":"(2025-03-22T00:35:23Z,2025-03-23T00:35:23Z]",
      "negative_certified":False,"source_first_public_availability_proven":False,
      "science_gate":"G3_NOT_ACCEPTED","original_files":[]}
 for sat,day,name in TARGETS:
  directory=f"{BASE}/goes{sat}/l2/data/xrsf-l2-avg1m_science/2025/03/"
  obj={"satellite":sat,"date_utc":f"2025-03-{day}","source_directory":directory,
       "initial_expected_name":name}
  try:
   index=get(directory,max_bytes=250_000).decode("utf8","replace")
   obj["source_directory_index_sha256"]=hashlib.sha256(index.encode("utf8")).hexdigest()
   matches=sorted(set(re.findall(r"sci_xrsf-l2-avg1m_g"+str(sat)+r"_d202503"+day+r"_v[0-9-]+\.nc",index)))
   if not matches: raise ValueError("NO_MATCHING_DATE_IN_LIVE_SOURCE_DIRECTORY")
   actual_name=matches[-1]
   url=directory+actual_name
   obj["original_url"]=url
   obj["actual_discovered_name"]=actual_name
   raw=get(url)
   if not(raw.startswith(b"CDF") or raw.startswith(b"\x89HDF")):
    raise ValueError("NOT_VALID_NETCDF_MAGIC")
   path=OUT/name
   path.write_bytes(raw)
   obj.update(state="BYTES_ACQUIRED",size_bytes=len(raw),
              sha256=hashlib.sha256(raw).hexdigest())
   try:obj["metadata"]=inspect(path);obj["state"]="BYTES_AND_METADATA_PASS"
   except Exception as e:obj["inspection_error"]=type(e).__name__+":"+str(e)[:300]
  except Exception as e:obj.update(state="SOURCE_NOT_RECOVERED",error=type(e).__name__+":"+str(e)[:350])
  out["original_files"].append(obj)
  (OUT/"original_xrs_nc_manifest.json").write_text(json.dumps(out,indent=2),encoding="utf8")
  print(json.dumps({k:obj.get(k) for k in ("satellite","date_utc","state","size_bytes","error")}),flush=True)
 print("RESULT",sum(x["state"]!="SOURCE_NOT_RECOVERED" for x in out["original_files"]),"/4",flush=True)
