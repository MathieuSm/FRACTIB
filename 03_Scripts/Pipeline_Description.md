# FRACTIB Pipeline

## HR-pQCT Scans
- Scan 5 stacks tibiae
- hFE standard evaluation of 4 stacks
- Transfer AIMs to hard drive (filezilla transfer type: binary)
    - AIMs to transfer: Measurement, SEG, CORT_MASK, TRAB_MASK
- Transfer AIMs to uCT (filezilla transfer type: binary)
- Run IPL script to uncompressed AIMs with common sizes
- Transfer AIMs to corresponding folder in 'FRACTIB/02_Data/01_HRpQCT/'
    - 1 folder per sample
    - Remove version extension (;1)

## Mechanical Tests
- Cut samples to similar region as hFE standard analyzed region
- Lap sample surface to get correct heigth, parallel and flat surfaces
- Scan sample with uCT (19.5um resolution)
- Test samples under compression up to failure, record max displacement -> record 3D displacements
- Scan tested sample with uCT (19.5um resolution) -> try to align position with first uct scan

## hFE analysis
- Register HRp-QCT on uCT scan using masks and store transformation parameters
- Crop masks uCT and HRp-QCT to have maximum common region with parallel surfaces
- Perform hFE using cropped masks. Set displacement to max experimental displacement
- Extract deformation gradient from simulation and decompose it to J and F_tilde
- Transform J and F_Tilde into uCT space

## Registration
- Register pre/post-test uct scans (rigid + bspline registrations) using cropped mask (common region with HRp-QCT)
- Extract deformation gradient, and decompose it to J and F_tilde

## High strain correlation
- Open J and F_tilde and for both HRp-QCT (hFE) and uCT file
- Perform linear regression with the 2 arrays of values (uCT vs hFE)
