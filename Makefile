IDIR =./include
CC=gcc
CFLAGS=-g -Wall -O2 -Wno-return-type -Wno-unused-variable -Wno-unused-function -I$(IDIR)
LFLAGS= -lz -lm 

ODIR=./obj/
SDIR=./src/
BDIR=./bin/

LIBS=-lm


#PREFIX is environment variable, but if it is not set, then set default value
ifeq ($(PREFIX),)
    PREFIX := /usr/local/bin
endif

install: $(BDIR)search-mgi-adapters $(BDIR)filter_reads_with_n
	install -d $(DESTDIR)$(PREFIX)
	install -m 755 $^ $(DESTDIR)$(PREFIX)

objects := $(ODIR)print-sequences.o \
	$(ODIR)store-primers.o $(ODIR)seqs_to_ints.o \
	$(ODIR)rob_dna.o $(ODIR)test.o $(ODIR)filter_reads_with_n.o \
	$(ODIR)match-paired-snps.o $(ODIR)search-paired-snp.o $(ODIR)trim-paired-snp.o

$(objects): $(ODIR)%.o: $(SDIR)%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@ $(FLAGS)


$(BDIR)test: $(ODIR)store-primers.o $(ODIR)seqs_to_ints.o $(ODIR)print-sequences.o $(ODIR)test.o $(ODIR)rob_dna.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

test: $(BDIR)test

$(BDIR)filter_reads_with_n: $(ODIR)rob_dna.o $(ODIR)filter_reads_with_n.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

filter_reads_with_n: $(BDIR)filter_reads_with_n


$(BDIR)trim-mgi-adapters: $(ODIR)seqs_to_ints.o $(ODIR)rob_dna.o $(ODIR)store-primers.o $(ODIR)match-paired-snps.o $(ODIR)trim-paired-snp.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

trim-mgi-adapters: $(BDIR)trim-mgi-adapters

$(BDIR)search-mgi-adapters: $(ODIR)seqs_to_ints.o $(ODIR)rob_dna.o $(ODIR)store-primers.o $(ODIR)match-paired-snps.o $(ODIR)search-paired-snp.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

search-mgi-adapters: $(BDIR)search-mgi-adapters

EXEC=filter_reads_with_n search-mgi-adapters trim-mgi-adapters test
all: $(addprefix $(BDIR), $(EXEC))



.PHONY: clean

clean:
	rm -fr bin/ obj/


