import os
import zipfile
import xml.etree.ElementTree as ET
import requests
from biom import load_table
import pandas as pd

def fetch_term_mapping():
    MDT_LOOKUP_URL = "https://mdt.gbif-test.org/service/terms/dwc"
    resp = requests.get(MDT_LOOKUP_URL)
    resp.raise_for_status()
    data = resp.json()
    mapping = {}
    for entry in data.values():
        qual = entry.get("qualName", "")
        name = entry.get("name", "")
        if qual and name:
            mapping[qual] = name
            if "/" in qual:
                short = qual.rsplit("/", 1)[-1]
                mapping[short] = name
    return mapping

def extract_term_value_pairs(xml_file: str, term_map: dict) -> list[tuple[str, str]]:
    ns = {"dwc": "http://rs.tdwg.org/dwc/text/"}
    tree = ET.parse(xml_file)
    root = tree.getroot()
    pairs = []
    for field in root.findall(".//dwc:field", ns):
        term_uri = field.attrib.get("term")
        default_value = field.attrib.get("default")
        if not term_uri or default_value is None:
            continue
        resolved_name = term_map.get(term_uri) or term_map.get(term_uri.rsplit("/", 1)[-1])
        if resolved_name:
            pairs.append((resolved_name, default_value))
        else:
            print(f"⚠️ Skipping unmapped term: {term_uri}")
    return pairs

def write_study_tsv(pairs: list[tuple[str, str]], filename="Study.tsv"):
    with open(filename, "w", encoding="utf-8") as f:
        f.write("term\tvalue\n")
        for term, value in pairs:
            f.write(f"{term}\t{value}\n")

def biom_to_tables(biom_file):
    table = load_table(biom_file)
    otu_df = table.to_dataframe(dense=True)
    sample_metadata = {
        sid: table.metadata(sid, axis='sample') for sid in table.ids(axis='sample')
        if table.metadata(sid, axis='sample')
    }
    sample_df = pd.DataFrame.from_dict(sample_metadata, orient='index')
    sample_df = sample_df.drop(columns=['further terms', 'readCount', 'id'], errors='ignore')
    sample_df = sample_df.reset_index().rename(columns={'index': 'id'})
    obs_metadata = {
        oid: table.metadata(oid, axis='observation') for oid in table.ids(axis='observation')
        if table.metadata(oid, axis='observation')
    }
    taxonomy_df = pd.DataFrame.from_dict(obs_metadata, orient='index')
    if 'taxonomy' in taxonomy_df.columns:
        taxonomy_df['taxonomy'] = taxonomy_df['taxonomy'].apply(
            lambda x: "; ".join(x) if isinstance(x, list) else x
        )
    taxonomy_df = taxonomy_df.drop(columns=['further terms', 'TaxonID', 'id'], errors='ignore')
    taxonomy_df = taxonomy_df.reset_index().rename(columns={'index': 'id'})
    return otu_df, sample_df, taxonomy_df

def main():
    biom_file = next((f for f in os.listdir() if f.endswith(".h5") or f == "data.biom.h5"), None)
    if not biom_file:
        print("❌ No BIOM file found (expected .h5 file or 'data.biom.h5').")
        return
    otu_df, sample_df, taxonomy_df = biom_to_tables(biom_file)
    otu_df.to_csv("OTU_table.csv")
    sample_df.to_csv("Samples.csv", index=False)
    taxonomy_df.to_csv("Taxonomy.csv", index=False)
    print("✅ BIOM tables written to OTU_table.csv, Samples.csv, and Taxonomy.csv")

    xml_path = None
    if os.path.isdir("archive"):
        xml_path = os.path.join("archive", "meta.xml")
    elif os.path.isfile("archive.zip"):
        with zipfile.ZipFile("archive.zip", "r") as zip_ref:
            zip_ref.extractall("archive")
        xml_path = os.path.join("archive", "meta.xml")
    if not xml_path or not os.path.isfile(xml_path):
        print("❌ meta.xml not found in archive or archive.zip.")
        return

    term_mapping = fetch_term_mapping()
    term_value_pairs = extract_term_value_pairs(xml_path, term_mapping)
    write_study_tsv(term_value_pairs)
    print(f"✅ Extracted {len(term_value_pairs)} term-value pairs to Study.tsv")

if __name__ == "__main__":
    main()
