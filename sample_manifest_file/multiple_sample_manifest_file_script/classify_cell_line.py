#!/usr/bin/env python
import sys
import re

for line in sys.stdin:
	lower_line = line.lower()
	notrna = 'FALSE'
	reference = 'FALSE'
	small = 'FALSE'
	differentiated = 'FALSE'
	if 'chipseq' in lower_line:
		notrna = 'TRUE'
	if 'mirna' in lower_line or 'srna' in lower_line \
		or 'small rna' in lower_line or 'microrna' in lower_line:
		small = 'TRUE'
	single = 'FALSE'
	if 'single-cell' in lower_line or 'single cell' in lower_line:
		single = 'TRUE'
	if '2011-11-22T10:02:05' in line or '2010-04-21T16:00' in line:
		small = 'TRUE'
	if 'Skin psoriasis involved' in line \
		or 'Skin psoriasis uninvolved' in line \
		or 'Skin Normal' in line:
		small = 'TRUE'
	if 'METSIM' in line:
		small = 'TRUE'
	if '2014-03-06T16:48:' in line \
		and 'source_name:Dorsolateral prefrontal cortex;' in line:
		small = 'TRUE'
	if '2014-03-06T16:48:' in line and ';subject id:' in line:
		small = 'TRUE'
	if single =='TRUE' or small == 'TRUE' or notrna == 'TRUE':
		#print '\t'.join([notrna, single, small, reference,
		#					differentiated, 'NA', 'NA', line]),
		continue
	cell_line = lower_line.strip().replace(
						'line:', '\x1c'
					).partition('\x1c')[2].partition(';')[0]
	# Cover exceptions
	cell_line = cell_line.strip()
	if cell_line == 'yes' and 'HCC827' in line:
		cell_line = 'HCC827'
	if cell_line == 'H1 hESCs':
		cell_line = 'H1'
	if cell_line.endswith('cells'):
		cell_line = cell_line[:-5].strip()
	if cell_line.endswith('cell'):
		cell_line = cell_line[:-4].strip()
	if cell_line.startswith('breast cancer cell line'):
		cell_line = cell_line[23:].strip()
	if cell_line in [str(number) for number in 
							[10847,12878,12890,12891,12892,18486,
								18505,18526,18951,19099,19193,19238,
								19239,19240,2255,2588,2610,2630]
						] and 'lymphoblastoid' in lower_line:
		cell_line = 'LCL'
	if cell_line.startswith('NA') or cell_line.startswith('GM'):
		try:
			int(cell_line[2:])
		except ValueError:
			pass
		else:
			if len(cell_line) == 7:
				cell_line = 'LCL'
	if 'geuv' in lower_line and 'lymphoblastoid' in lower_line:
		cell_line = 'LCL'
	if cell_line == '1221' and 'NUT' in line:
		'''patient-derived cell line: http://www.ncbi.nlm.nih.gov/
		geo/query/acc.cgi?acc=GSE57222'''
		cell_line = 'patient-derived; lung metastasis'
	if cell_line in ['141', '510'] and 'liposarcoma' in line:
		cell_line = 'LPS' + cell_line
	if cell_line == 'A172, A431, D538MG':
		cell_line = '+'.join([
						a_line.strip() for a_line in cell_line.split(',')
					])
	if 'HGDP' in cell_line and 'lymphoblast' in lower_line:
		cell_line = 'LCL'
	cell_line = cell_line.strip()
	if 't24' in lower_line and 'bladder carcinoma cell' in lower_line:
		cell_line = 'T24'
	if not cell_line and ('hela cell' in lower_line
							or '_hela' in lower_line
							or 'hela_' in lower_line
							or ':hela' in lower_line
							or 'helamock' in lower_line
							or 'srs259977' in lower_line):
		cell_line = 'HeLa'
	if not cell_line and 'lung adenocarcinoma cell lines' in lower_line:
		cell_line = line.replace(
							'lung adenocarcinoma cell lines', '\x1c'
						).partition('\x1c')[0].rpartition('\t')[2].strip()
	if not cell_line and 'ovarian cancer cell line' in lower_line:
		cell_line = line.replace(
							'ovarian cancer cell line', '\x1c'
						).partition('\x1c')[0].rpartition(';')[2].strip()
	if not cell_line and 'Breast Cancer Cell Line' in line:
		cell_line = line.split('\t')[4]
	if not cell_line and 'lymphoblastoid cell line' in lower_line:
		cell_line = 'LCL'
	if not cell_line and ('Leukaemic cell line treated for 24h with vehicle'
							in line):
		cell_line = '-'.join(line.split('\t')[5].split('_')[:2])
	if not cell_line and 'prostate cancer cell line' in line:
		cell_line = line.split('\t')[4].split(' ')[0]
	if not cell_line and 'T24_' in line and 'bladde cancer cell' in line:
		cell_line = 'T24'
	if not cell_line and 'A549' in line:
		cell_line = 'A549'
	if not cell_line and 'hesc' in lower_line:
		cell_line = 'hESC'
	if 'glioma tissue -' in lower_line:
		cell_line = 'FALSE'
	if not cell_line and 'hTERT-RPE1' in line:
		cell_line = 'hTERT-RPE1'
	if not cell_line and ('This sample is the RNA that isolated from '
		                  'peripheral blood samples') in line:
		cell_line = 'FALSE'
	if not cell_line and 'UHRR' in line and 'HBRR' in line:
		cell_line = 'FALSE'
		reference = 'UHRR+HBRR'
	if not cell_line and 'kasumi1' in lower_line:
		cell_line = 'Kasumi-1'
	if not cell_line and 'primary cell' in lower_line:
		cell_line = 'FALSE'
	if not cell_line and 'Total RNAs of a human adult brain' in line:
		cell_line = 'FALSE'
	if not cell_line and '(contact:Al Forrest)' in line:
		cell_line = 'FALSE'
	if not cell_line and 'source_name:Non-smoker female' in line:
		cell_line = 'FALSE'
	if not cell_line \
		and 'source_name:Breast tissue;tissue:Breast;disease status' in line:
		cell_line = 'FALSE'
	if not cell_line and 'VADS_BW' in line:
		cell_line = 'FALSE'
	if not cell_line and 'E-MTAB-513' in line:
		cell_line = 'FALSE'
	if not cell_line and '(HBRR)' in line:
		cell_line = 'FALSE'
		reference = 'HBRR'
	if not cell_line and 'collection site:Atlanta VA Medical Center' in line:
		# prostate
		cell_line = 'FALSE'
	if not cell_line and 'whole_blood_derived_from_typhoid_challenge' in line:
		cell_line = 'FALSE'
	if not cell_line and 'lncap' in lower_line:
		cell_line ='LNCaP'
	if not cell_line and 'Huh7 human hepatoma cells' in line:
		cell_line = 'Huh7'
	if not cell_line and 'pancreatic islets' in line:
		cell_line = 'FALSE'
	if not cell_line and '(UHRR)' in line:
		cell_line = 'FALSE'
		reference = 'UHRR'
	if not cell_line and line.split('\t')[1].startswith('2012-01-10T13:13:0'):
		cell_line = 'FALSE'
	if not cell_line and 'ASE_hiseq' in line or 'ASE_miseq' in line:
		cell_line = 'FALSE'
	if not cell_line and 'source_name:renal;' in line:
		cell_line = 'FALSE'
	if not cell_line and 'specimen with known storage state' in line:
		cell_line = 'FALSE'
	if not cell_line and 'A431' in line:
		cell_line = 'A431'
	if not cell_line and 'U-251' in line:
		cell_line = 'U-251'
	if not cell_line and 'U2OS' in line:
		cell_line = 'U2OS'
	if not cell_line and ('A-549' in line or 'A549' in line):
		cell_line = 'A-549'
	if not cell_line \
		and 'source_name:Neural differentiation;' in line:
		cell_line = 'H1'
		differentiated = line.split('\t')[4].partition('_')[0]
	if not cell_line and 'human oral biofilm' in line:
		cell_line = 'FALSE'
	if not cell_line and 'Large intestine tumor' in line:
		cell_line = 'FALSE'
	if not cell_line and 'source_name:non-neoplastic brain tissue' in line:
		cell_line = 'FALSE'
	if not cell_line and 'source_name:skeletal muscle;' in line:
		cell_line = 'FALSE'
	if not cell_line and 'GSM754335' in line:
		cell_line = 'LCL'
	if not cell_line and 'glioblastoma cell-line' in line:
		cell_line = line.split('\t')[5].split(' ')[0]
	if not cell_line and ';biomaterial_provider:Nijman' in line:
		cell_line = line.split('\t')[4][4:].replace('.', ' ')
	if not cell_line and 'HBRR_' in line:
		cell_line = 'FALSE'
		reference = 'HBRR'
	if not cell_line and ' cells' in lower_line:
		if 'peripheral blood mononuclear ' in lower_line \
			or 'source tissue:tonsil' in lower_line \
			or 'myeloid cells' in lower_line \
			or 'CD34+' in line or 'CD4' in line or 'Treg' in line \
			or 'source_name:Large airway epithelial cells;' in line \
			or 'source_name:neonatal dermal BECs' in line \
			or 'Induced pluripotent' in line \
			or 'source_name:whole peripheral blood;' in line \
			or 'cell type:Primary Dermal Fibroblast;' in line \
			or 'airway basal' in line:
			cell_line = 'FALSE'
		else:
			pass
	if not cell_line:
		print line,