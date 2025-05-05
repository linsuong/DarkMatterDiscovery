import re
import pandas as pd

with open('exlcusion_limits.tex', 'r') as f:
    content = f.read()

# Improved pattern to capture all number formats and scientific notation
row_pattern = re.compile(
    r'(\d+)\s*&\s*(\d+)\s*&\s*(\d+)\s*&\s*([\d\.]+(?:\\times10\^\{\d+\})?|-)\s*&\s*(\d+|-)(?:[^\\]|$)', 
    re.MULTILINE
)

data = []
for row in row_pattern.finditer(content):
    m_D1 = int(row.group(1))
    delta_m_plus = int(row.group(2))
    delta_m_0 = int(row.group(3))
    
    # Process cross-section
    sigma = row.group(4)
    if sigma == '-':
        continue  # Skip rows with missing cross-sections
    if '\\times10' in sigma:
        sigma = sigma.replace('\\times10^{', 'e').replace('}', '')
    
    # Process error
    err = row.group(5)
    if err == '-':
        continue  # Skip rows with missing errors
    
    data.append([
        m_D1, delta_m_plus, delta_m_0,
        float(sigma), float(err)
    ])

# Create DataFrame
df = pd.DataFrame(data, columns=[
    'm_D1', 'delta_m_plus', 'delta_m_0',
    'sigma_2l_A', 'err_2l_A'
])

print(f"Total extracted points: {len(df)}")
df.to_csv('cross_sections/limits_extracted_full.dat', sep='\t', index=False, float_format='%.6g')