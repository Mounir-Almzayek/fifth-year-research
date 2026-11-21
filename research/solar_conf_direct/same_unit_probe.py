#!/usr/bin/env python3
"""Alternative read-only retrieval of conf_disambig by original JSOC SUM folder.
The source manifest establishes the source record folder using BR/BITMAP paths;
the conf segment suffix is a discovery candidate, NOT independently DRMS-confirmed.
"""
import hashlib, io, json, pathlib, unittest
import requests,numpy as np
from astropy.io import fits
S=requests.Session();S.headers["User-Agent"]="Solar-Intelligence-research-source-recovery/1.0"
OUT=pathlib.Path("solar-conf-direct-evidence");OUT.mkdir(exist_ok=True)
RECORDS=[
(7115,"/SUM98/D978177825/S00000",(448,688),"2017.09.04_00:00:00_TAI"),
(7122,"/SUM95/D978180834/S00000",(317,549),"2017.09.04_00:00:00_TAI"),
(9781,"/SUM82/D1666182645/S00000",(461,1061),"2023.07.16_00:36:00_TAI"),
(9802,"/SUM7/D1666183096/S00000",(732,1947),"2023.07.16_00:36:00_TAI"),
(11149,"/SUM97/D1753927374/S00000",(513,1202),"2024.05.10_00:00:00_TAI"),
(11170,"/SUM88/D1758236239/S00000",(415,771),"2024.05.13_00:36:00_TAI"),
(11220,"/SUM89/D1758228691/S00000",(96,132),"2024.05.13_00:36:00_TAI"),
(12941,"/SUM87/D1869862244/S00000",(392,1031),"2025.03.22_00:36:00_TAI")]
def probe(h,folder,shape,epoch):
 rec={"harp":h,"original_record_t_rec":epoch,"source_folder_from_original_manifest":folder,
     "segment_path_DRMS_confirmed":False,"status":"NOT_RECOVERED","attempts":[]}
 for host in ("https://jsoc1.stanford.edu","https://jsoc.stanford.edu"):
  url=host+folder+"/conf_disambig.fits"
  try:
   with S.get(url,timeout=(8,17),stream=True) as resp:
    resp.raise_for_status()
    content=bytearray()
    for buf in resp.iter_content(131072):
     content.extend(buf)
     if len(content)>12_000_000:raise ValueError("FILE_SIZE_LIMIT")
   raw=bytes(content)
   if not raw.startswith(b"SIMPLE  =") or len(raw)%2880:raise ValueError("INVALID_FITS_ENVELOPE")
   (OUT/f"harp_{h}_conf_disambig.fits").write_bytes(raw)
   rec.update(status="FITS_BYTES_RECOVERED_QA_PENDING",source_url=url,
              sha256=hashlib.sha256(raw).hexdigest(),size_bytes=len(raw),
              nominal_source_shape=list(shape))
   try:
    with fits.open(io.BytesIO(raw),checksum=True,memmap=False) as hdus:
     imgs=[hdu for hdu in hdus if getattr(hdu.data,"ndim",0)==2]
     if len(imgs)!=1:raise ValueError("EXPECTED_ONE_FITS_IMAGE")
     image=np.asarray(imgs[0].data)
     if tuple(image.shape)!=tuple(shape):raise ValueError("GRID_SHAPE_MISMATCH")
     c=[hd.verify_checksum() for hd in hdus if "CHECKSUM" in hd.header]
     d=[hd.verify_datasum() for hd in hdus if "DATASUM" in hd.header]
     if not c or not d or any(k!=1 for k in c+d):raise ValueError("BAD_FITS_CHECKSUM")
     finite=image[np.isfinite(image)]
     vals,counts=np.unique(finite,return_counts=True)
     values={str(int(i)):int(n) for i,n in zip(vals,counts)}
     rec.update(status="ORIGINAL_FITS_QA_PASS",shape=list(image.shape),
                unique_values=values,n90=int(np.count_nonzero(image==90)),
                invalid_pixel_count=int(image.size-finite.size),
                checksum=c,datasum=d)
   except Exception as qa_err:
    rec["qa_error"]=type(qa_err).__name__+":"+str(qa_err)[:300]
   break
  except Exception as exc:
   rec["attempts"].append({"url":url,"error":type(exc).__name__+":"+str(exc)[:160]})
 return rec
class Checks(unittest.TestCase):
 def test_eight_unique_sources(self):self.assertEqual(len({a[0] for a in RECORDS}),8)
 def test_eight_storage_units(self):self.assertEqual(len({a[1] for a in RECORDS}),8)
 def test_target_sample(self):self.assertEqual(RECORDS[0][0],7115)
if __name__=="__main__":
 res=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Checks))
 if not res.wasSuccessful():raise SystemExit(1)
 rows=[]
 for args in RECORDS:
  rec=probe(*args);rows.append(rec)
  (OUT/"same_unit_direct_probe_results.json").write_text(json.dumps({"method":"SAME_SUM_STORAGE_UNIT_WITH_CANDIDATE_CONF_SEGMENT_PATH","rows":rows,"G1":"OPEN"},indent=2),encoding="utf-8")
  print(rec["harp"],rec["status"],rec.get("sha256",""),flush=True)
 print("RECOVERED",sum(x["status"]!="NOT_RECOVERED" for x in rows),"/8",flush=True)
