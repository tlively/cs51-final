# This is the Makefile for building the project! yay!

all: game
	@echo "building all!"

game: game.c physics.o graphics.o
	gcc -o game game.c -lm

physics.o: physics.c
	gcc -o physics.o -c physics.c -lm

graphics.o: graphics.c
	gcc -o graphics.o -c graphics.c `sld-config --cflags --libs` -lm

clean:
	rm -f game physics.o graphics.o

.PHONY: all
