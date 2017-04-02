.PHONY: all

all: readme.txt
	@true

readme.txt:
	@echo "DOWNLOADING DATASET"
	@wget -c http://imagecompression.info/test_images/gray8bit.zip
	@unzip gray8bit.zip
	@rm gray8bit.zip
	#for entropy encoding, it is more representative to focus on real images
	@rm zone_plate.pgm artificial.pgm

clean:
	rm $(wildcard *.pgm) readme.txt
