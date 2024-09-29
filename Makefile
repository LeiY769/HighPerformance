CC=gcc
CFLAGS=
LDFLAGS= -lm -O3
LD=gcc

NAME = main
MODULES = main.c \
		  utils.c

OBJECTS = $(MODULES:.c=.o)

all: $(NAME)

$(NAME) : $(OBJECTS)
	$(LD) $(OBJECTS) -o $(NAME) $(LDFLAGS)

$(OBJECTS): %.o: %.c 
	$(CC) -c $? -o $@ $(CFLAGS)

clean:
	rm -rf *.o *.vti *.pvd $(NAME)

run: $(NAME)
	./$(NAME) param_simple.txt

.PHONY: all clean run