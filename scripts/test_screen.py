import requests

csv = """gene,guide,control_count,treatment_count
ESS001,ESS001_g1,1000,50
ESS001,ESS001_g2,800,40
ESS001,ESS001_g3,900,45
NEU001,NEU001_g1,500,490
NEU001,NEU001_g2,600,610
NEU001,NEU001_g3,550,545
"""

resp = requests.post(
    'http://localhost:8000/api/screen/analyze',
    files={'file': ('test.csv', csv.encode(), 'text/csv')},
    data={
        'guide_col': 'guide',
        'gene_col': 'gene',
        'control_col': 'control_count',
        'treatment_col': 'treatment_count',
        'min_reads': '30'
    }
)
print('Status:', resp.status_code)
print('Response:', resp.text)
