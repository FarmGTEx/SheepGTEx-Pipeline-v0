#!/bin/bash
source activate HYFA

cd ../..
python train_sheepgoat_gtex.py --config configs/default.yaml \
	--expressfile sheep/v4.min25.split/express.csv \
	--metafile sheep/v4.min25.split/meta.txt \
	--tiscolorfile sheep/tissue_color.list \
	--mintissuenum 2 --minsamplenum 25 \
	--traindonorfile sheep/v4.min25.split/splits/train_tis2.txt \
	--valdonorfile sheep/v4.min25.split/splits/val_tis2.txt \
	--testdonorfile sheep/v4.min25.split/splits/test_tis2.txt \
	--modelfile sheep/v4.min25.split/sheep_tis2.model.pth

