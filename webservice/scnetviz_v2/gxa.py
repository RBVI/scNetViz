import requests, zipfile, io, os

GXA_MTX_URI = "https://www.ebi.ac.uk/gxa/sc/experiment/%s/download/zip?fileType=normalised"
FILE = "%s.aggregated_filtered_normalised_counts"

def fetch_experiment(accession, cachePath):
  url = GXA_MTX_URI%accession
  r = requests.get(url)
  if r.ok:
    z = zipfile.ZipFile(io.BytesIO(r.content))
    z.extractall(cachePath)

    # Now, rename the files
    file_prefix = cachePath+'/'+FILE%accession
    os.rename(file_prefix+".mtx", cachePath+"/matrix.mtx")
    os.rename(file_prefix+".mtx_cols", cachePath+"/barcodes.tsv")
    os.rename(file_prefix+".mtx_rows", cachePath+"/genes.tsv")

  return r.status_code
    
