#! /usr/bin python

from hFE_6ReadResults import Main as hFE_6ReadResults

Samples = ['432_L_77_F',
           '433_R_77_F',
           '434_L_90_F',
           '435_R_90_F',
           '436_L_90_F',
           '437_R_90_F',
           '438_L_71_F',
           '439_R_71_F',
           '440_L_64_M',
           '441_R_64_M',
           '442_R_75_F',
           '443_L_73_F',
           '444_R_92_F',
           '445_R_93_F',
           '446_R_75_F',
           '447_L_83_M',
           '448_L_80_M',
           '449_L_93_F',
           '450_L_77_F',
           '451_L_75_F',
           '452_L_75_F',
           '453_R_79_M',
           '454_L_94_F',
           '455_L_97_F',
           '456_R_97_F']

class Arguments():

	def __init__(self):
		self.Folder = 'FRACTIB'

Arguments = Arguments()

for Sample in Samples:

	Arguments.Sample = Sample

	hFE_6ReadResults(Arguments)

