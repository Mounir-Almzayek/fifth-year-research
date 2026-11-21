#!/usr/bin/env python3
"""Read-only G1 source probe: 8 exact definitive CEA conf_disambig FITS.
No gate acceptance, no prediction, no credentials; preserve source bytes and failures.
"""
import csv, hashlib, io, json, pathlib, re, time, unittest
from urllib.parse import urlparse
import requests
import numpy as np
from astropy.io import fits

SERIES="hmi.sharp_cea_720s"
RECORDS=[
    (7115,"2017.09.04_00:00:00_TAI","12673"),
    (7122,"2017.09.04_00:00:00_TAI","12677"),
    (9781,"2023.07.16_00:36:00_TAI","13363"),
    (9802,"2023.07.16_00:36:00_TAI","13372"),
    (11149,"2024.05.10_00:00:00_TAI","13664,13668"),
    (11170,"2024.05.13_00:36:00_TAI","13667"),
    (11220,"2024.05.13_00:36:00_TAI","13677"),
    (12941,"2025.03.22_00:36:00_TAI","14035"),
]
HOSTS=("https://jsoc.stanford.edu","https://jsoc1.stanford.edu","http://jsoc.stanford.edu")
KEYS="HARPNUM,T_REC,T_OBS,NOAA_ARS,NOAA_NUM,QUALITY,CMASK,USFLUX,AREA_ACR,CDELT1,CDELT2"
OUT=pathlib.Path("solar-conf-evidence")
OUT.mkdir(parents=True,exist_ok=True)
S=requests.Session()
S.headers["User-Agent"]="Solar-Intelligence-G1-Original-Conf-Read-Only/1.0"

def get(url,params=None,limit=8_000_000):
    with S.get(url,params=params,timeout=(10,35),stream=True) as response:
        response.raise_for_status()
        data=bytearray()
        for part in response.iter_content(131072):
            data.extend(part)
            if len(data)>limit: raise ValueError("SOURCE_SIZE_LIMIT")
        return bytes(data),response.url
def vals(columns,count):
    result=[{} for _ in range(count)]
    for c in columns:
        if len(c.get("values",[]))!=count: raise ValueError("BAD_COLUMN_LENGTH")
        for i,v in enumerate(c["values"]):result[i][c["name"]]=v
    return result
def one_record(doc):
    if doc.get("status")!=0 or doc.get("count")!=1:raise ValueError("INVALID_SOURCE_ROW_COUNT_OR_STATUS")
    return vals(doc.get("keywords",[]),1)[0],vals(doc.get("segments",[]),1)[0]
def safe_source(path):
    return bool(isinstance(path,str) and re.fullmatch(r"/SUM[0-9]+/D[0-9]+/S[0-9]+/conf_disambig\.fits",path))
def source_query(harp,epoch):
    params={"op":"rs_list","ds":f"{SERIES}[{harp}][{epoch}]","key":KEYS,"seg":"conf_disambig"}
    attempts=[]
    for host in HOSTS:
        try:
            raw,url=get(host+"/cgi-bin/ajax/jsoc_info",params,limit=2_000_000)
            doc=json.loads(raw)
            row,segments=one_record(doc)
            if row.get("HARPNUM")!=str(harp) or row.get("T_REC")!=epoch:raise ValueError("SOURCE_IDENTITY_MISMATCH")
            return {"ok":True,"url":url,"raw":raw,"row":row,"segment":segments.get("conf_disambig"),"attempts":attempts}
        except Exception as ex:
            attempts.append({"url":host,"error":type(ex).__name__+":"+str(ex)[:140]})
    return {"ok":False,"attempts":attempts}
def fetch_fits(path,harp):
    attempts=[]
    if not safe_source(path):return {"ok":False,"reason":"NO_CONF_SEGMENT_OR_UNSAFE_PATH","path":path}
    for host in ("https://jsoc1.stanford.edu","https://jsoc.stanford.edu"):
        try:
            raw,url=get(host+path,limit=10_000_000)
            if len(raw)%2880 or not raw.startswith(b"SIMPLE  ="):
                raise ValueError("INVALID_FITS_ENVELOPE")
            with fits.open(io.BytesIO(raw),memmap=False,checksum=True) as hdus:
                images=[h for h in hdus if getattr(h.data,"ndim",0)==2]
                if len(images)!=1:raise ValueError("EXPECTED_ONE_IMAGE")
                obj=images[0]
                data=np.asarray(obj.data)
                if not np.isfinite(data).all():raise ValueError("NONFINITE_PIXEL")
                check=[h.verify_checksum() for h in hdus if "CHECKSUM" in h.header]
                sums=[h.verify_datasum() for h in hdus if "DATASUM" in h.header]
                if not check or not sums or any(v!=1 for v in check+sums):
                    raise ValueError("CHECKSUM_OR_DATASUM_FAILURE")
                hist={str(int(v)):int(c) for v,c in zip(*np.unique(data,return_counts=True))}
                result=dict(ok=True,sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw),
                            source_url=url,shape=list(data.shape),unique_values=hist,
                            checksum_verified=check,datasum_verified=sums,
                            n_high_confidence_90=int(np.count_nonzero(data==90)))
            (OUT/f"harp_{harp}_conf_disambig.fits").write_bytes(raw)
            return result
        except Exception as ex:
            attempts.append({"host":host,"error":type(ex).__name__+":"+str(ex)[:160]})
    return {"ok":False,"reason":"SEGMENT_HTTP_OR_FITS_FAILED","attempts":attempts}

class GuardTests(unittest.TestCase):
    def test_reject_hostile_path(self):
        self.assertFalse(safe_source("/etc/passwd"))
    def test_reject_other_segment(self):
        self.assertFalse(safe_source("/SUM12/D10/S00001/Br.fits"))
    def test_accept_source_path(self):
        self.assertTrue(safe_source("/SUM12/D10/S00001/conf_disambig.fits"))
    def test_not_accept_missing_records(self):
        with self.assertRaises(ValueError):one_record({"status":0,"count":0})
    def test_single_record_column(self):
        self.assertEqual(vals([{"name":"HARPNUM","values":["7115"]}],1)[0]["HARPNUM"],"7115")
    def test_fail_column_length(self):
        with self.assertRaises(ValueError):vals([{"name":"X","values":[]}],1)
if __name__=="__main__":
    tests=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(GuardTests))
    if not tests.wasSuccessful():raise SystemExit("TEST_GUARD_FAIL")
    output={"scientific_acceptance":"G1_OPEN","source_series":SERIES,
            "public_first_availability":"UNKNOWN","results":[]}
    for harp,epoch,noaa in RECORDS:
        rec={"harp":harp,"epoch":epoch,"retrospective_noaa_expected":noaa,
             "source_series":SERIES,"accepted_for_training":False}
        queried=source_query(harp,epoch)
        rec["query_attempts"]=queried["attempts"]
        if queried["ok"]:
            raw=queried.pop("raw")
            p=OUT/f"harp_{harp}_exact_source_query.json"
            p.write_bytes(raw)
            rec["source_query_sha256"]=hashlib.sha256(raw).hexdigest()
            rec["source_query_url"]=queried["url"]
            rec["keywords"]=queried["row"]
            rec["segment_path"]=queried["segment"]
            rec["fits"]=fetch_fits(queried["segment"],harp)
        else:rec["fits"]={"ok":False,"reason":"DRMS_QUERY_FAILED"}
        output["results"].append(rec)
        (OUT/"source_probe_results.json").write_text(json.dumps(output,indent=2),encoding="utf-8")
        print(json.dumps({"harp":harp,"query_ok":queried["ok"],
                          "fits_ok":rec["fits"]["ok"],
                          "sha256":rec["fits"].get("sha256")},sort_keys=True),flush=True)
    with (OUT/"source_probe_summary.csv").open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["harp","epoch","noaa","source_query_hash","CMASK","USFLUX","conf_sha256","bytes","shape","n90","state"])
        for a in output["results"]:
            k=a.get("keywords",{});f=a["fits"]
            w.writerow([a["harp"],a["epoch"],a["retrospective_noaa_expected"],a.get("source_query_sha256",""),
                        k.get("CMASK",""),k.get("USFLUX",""),f.get("sha256",""),f.get("bytes",""),
                        f.get("shape",""),f.get("n_high_confidence_90",""),
                        "SOURCE_FITS_CHECKED" if f["ok"] else "NOT_RECOVERED"])
    print("SUMMARY",sum(x["fits"]["ok"] for x in output["results"]),"/",len(RECORDS),"original FITS",flush=True)
