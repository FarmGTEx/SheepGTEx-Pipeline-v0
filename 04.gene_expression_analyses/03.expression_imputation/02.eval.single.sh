#!/bin/bash
sourcetissue=$1
source activate HYFA

cd ../..
python eval_sheepgoat_gtex.single.py --config configs/default.yaml \
	--modelfile sheep/v4.min25.split/sheep_tis2.model.pth \
	--expressfile sheep/v4.min25.split/express.csv \
	--metafile sheep/v4.min25.split/meta.txt \
	--tiscolorfile sheep/tissue_color.list \
	--mintissuenum 2 --minsamplenum 25 --sourcetissue ${sourcetissue} \
	--traindonorfile sheep/v4.min25.split/splits/train_tis2.txt \
	--valdonorfile sheep/v4.min25.split/splits/val_tis2.txt \
	--testdonorfile sheep/v4.min25.split/splits/test_tis2.txt \
	--resultsfile sheep/v4.min25.split/sheep_tis2.${sourcetissue}.csv
