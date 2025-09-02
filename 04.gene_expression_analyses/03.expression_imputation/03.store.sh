#!/bin/bash
source activate HYFA

cd ../..
python store_sheepgoat_gtex.py --config configs/default.yaml \
	--modelfile sheep/v4.min25.split/sheep_tis2.model.pth \
	--expressfile sheep/v4.min25.split/express.csv \
	--metafile sheep/v4.min25.split/meta.txt \
	--tiscolorfile sheep/tissue_color.list \
	--mintissuenum 2 --minsamplenum 25 \
	--observedfile sheep/v4.min25.split/sheep_observed_normalised.csv \
	--imputedfile sheep/v4.min25.split/sheep_imputed_normalised.csv 
