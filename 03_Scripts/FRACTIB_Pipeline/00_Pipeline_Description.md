FRACTIB Pipeline

Steps to perform:

HR-pQCT
1. Scan 5 stacks tibiae
1. hFE standard evaluation of 4 stacks
2. Transfer AIMs to hard drive (filezilla transfer type: binary)
    AIMs to transfer: Measurement, SEG, CORT_MASK, TRAB_MASK
    
IPL
3. Transfer AIMs to uCT (filezilla transfer type: binary)
4. Run IPL script to uncompressed AIMs with common sizes
5. Transfer AIMs to corresponding folder in 'FRACTIB/02_Data/02_HRpQCT/' (1 folder per sample)

hFE
6. Run AIM2MHD.py to remove version extension and write MHDs files
7. Register HRp-QCT on uCT scan using gray values and masks files, store transformation parameters
8. Crop HR-pQCT masks using uCT mask. Idea is to keep only common regions
9. Perform hFE using cropped masks. Set displacement to max experimental displacement
10. Extract deformation gradient from simulation and decompose it to J and F_tilde
11. Use registration transformation parameters to tranform J and F_tilde into uCT space

Sample testing
12. Cut samples to similar region as hFE standard analyzed region
13. Lap sample surface to get correct heigth, parallel and flat surfaces
14. Scan sample with uCT (19.5um resolution)
15. Test samples under compression up to failure, record max displacement -> record 3D displacements
16. Scan failed sample with uCT (19.5um resolution) -> try to align position with first uct scan
17. Register pre/post-test sample uct scans (rigid + bspline registrations)
18. Extract deformation gradient, and decompose it to J and F_tilde

High strain correlation
19. Open HR-pQCT J or F_tilde and corresponding uCT file
20. Multiply both image with uCT mask
21. Perform linear regression with the 2 arrays of values (uCT vs hFE)
