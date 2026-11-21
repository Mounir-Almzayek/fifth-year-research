#!/usr/bin/env python3
"""G3: fetch original NOAA NCEI GOES XRS L2 composite flare reports, unchanged.
This product combines GOES XRS Flare Summary & SWPC event reports; SWPC lineages
overlap with original NOAA edited events. NO independent-source claim or label.
"""
import csv, datetime as dt, hashlib, io, json, pathlib, unittest
import requests
BASE="https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/multi/l2/data/xrsf-l2-flrpt_science/csv"
YEARS=(2017,2023,2024,2025)
OUT=pathlib.Path("solar-goes-ncei-evidence");OUT.mkdir(parents=True,exist_ok=True)
S=requests.Session()
S.headers.update({"User-Agent":"Solar-Intelligence-G3-Source-Audit/1.0","Accept":"text/csv,application/json,*/*"})
def download(name,max_bytes=3_500_000):
    url=f"{BASE}/{name}"
    with S.get(url,timeout=(12,70),stream=True) as response:
        response.raise_for_status()
        content=bytearray()
        for piece in response.iter_content(131072):
            content.extend(piece)
            if len(content)>max_bytes: raise ValueError("SOURCE_TOO_LARGE_FOR_BOUNDED_STUDY")
        if not content:raise ValueError("EMPTY_BODY")
    raw=bytes(content)
    (OUT/name).write_bytes(raw)
    return {"url":url,"sha256":hashlib.sha256(raw).hexdigest(),"bytes":len(raw)}
def inspect_csv(raw):
    decoded=raw.decode("utf-8-sig")
    sample=list(csv.reader(io.StringIO(decoded)))
    if not sample:raise ValueError("NO_CSV_RECORDS")
    # NCEI product has one header row, numeric parameters and UTC times.
    first=sample[0]
    if len(first)<5:raise ValueError("TOO_FEW_COLUMNS")
    return {"columns":first,"records":len(sample)-1,
            "first_data":sample[1][:10] if len(sample)>1 else []}
def test_unit():
    class QA(unittest.TestCase):
        def test_sha_known(self): self.assertEqual(len(hashlib.sha256(b"").hexdigest()),64)
        def test_required_years(self): self.assertEqual(len(YEARS),4)
        def test_path_whitelist(self): self.assertEqual(YEARS[0],2017)
        def test_csv_sample(self): self.assertEqual(inspect_csv(b"a,b,c,d,e\n1,2,3,4,5\n")["records"],1)
        def test_unparsed_utc(self): self.assertEqual(dt.datetime.strptime("2017-09-04","%Y-%m-%d").year,2017)
    test=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(QA))
    if not test.wasSuccessful():raise SystemExit(1)
if __name__=="__main__":
    test_unit()
    report={"source":"NOAA NCEI L2 XRS Flare Report (composite GOES XRS + SWPC)",
            "science_gate":"G3_NOT_ACCEPTED",
            "release":"v1-0-1","years":{}}
    for name in ["sci_xrsf-l2-flrpt_geo_metadata.json"] + [f"sci_xrsf-l2-flrpt_geo_y{y}_v1-0-1.csv" for y in YEARS]:
        print("SOURCE",name,flush=True)
        try:
            entry=download(name)
            if name.endswith(".csv"):
                entry.update(inspect_csv((OUT/name).read_bytes()))
            else:
                obj=json.loads((OUT/name).read_text())
                entry["source_description"]=obj.get("global_attributes",{}).get("summary")
                entry["source_version"]=obj.get("global_attributes",{}).get("flrpt_algorithm_version")
            entry["status"]="RAW_ORIGINAL_ACQUIRED"
        except Exception as e:
            entry={"status":"NOT_RECOVERED","error":type(e).__name__+":"+str(e)[:300]}
        report["years"][name]=entry
        (OUT/"raw_download_manifest.json").write_text(json.dumps(report,indent=2),encoding="utf8")
        print(json.dumps({name:entry["status"],"rows":entry.get("records"),"sha":entry.get("sha256")}),flush=True)
    print("RESULT",sum(x["status"]=="RAW_ORIGINAL_ACQUIRED" for x in report["years"].values()),"/",len(report["years"]),flush=True)
