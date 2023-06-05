CFLAGS=-Wall 
LIBS=-O2 -L/usr/X11R6/lib -lm -lpthread -lX11
EXE=TASP

all:
	mkdir -p res
	g++ $(CFLAGS) TASP.cpp -o $(EXE) $(LIBS)

test:
	./TASP -i ./data/test_img.jpg -k 450 -m 0.1 -outm test_img_labelmap.png -outb test_img_border.png

test_list:
	./scripts/test_list.sh ./data/list_file.txt ./data/ 450 0.1
