IDIR =./include
CC=gcc
# note that usually use -O2 but for valgrind debugging use -O0 which is slower but more accurate
CFLAGS=-g -Wall -O0 -Wno-return-type -Wno-unused-variable -Wno-unused-function -I$(IDIR)
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

OBJ=print-sequences store-primers seqs_to_ints rob_dna test filter_reads_with_n match-paired-snps search-paired-snp trim-paired-snp \
	create-snps hash read_primers match-all-snps search-adapter-file primer-match-counts
objects := $(addsuffix .o, $(addprefix $(ODIR), $(OBJ)))

$(objects): $(ODIR)%.o: $(SDIR)%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@ $(FLAGS)

BASE=create-snps hash rob_dna seqs_to_ints store-primers 

TEST=${BASE} print-sequences test
testobj := $(addsuffix .o, $(addprefix $(ODIR), $(TEST)))
$(BDIR)test: $(testobj)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

test: $(BDIR)test


FRN=rob_dna filter_reads_with_n
frnobj := $(addsuffix .o, $(addprefix $(ODIR), $(FRN)))
$(BDIR)filter_reads_with_n: $(frnobj)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)
filter_reads_with_n: $(BDIR)filter_reads_with_n


TRIM=${BASE} match-paired-snps trim-paired-snp
trimobj := $(addsuffix .o, $(addprefix $(ODIR), $(TRIM)))
$(BDIR)trim-mgi-adapters: $(trimobj)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)
trim-mgi-adapters: $(BDIR)trim-mgi-adapters


SEARCHMGI=${BASE} match-paired-snps search-paired-snp 
searchmgiobj := $(addsuffix .o, $(addprefix $(ODIR), $(SEARCHMGI)))
$(BDIR)search-mgi-adapters: $(searchmgiobj)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)
search-mgi-adapters: $(BDIR)search-mgi-adapters


SEARCHAF=${BASE} match-all-snps primer-match-counts read_primers search-adapter-file 
searchafobj := $(addsuffix .o, $(addprefix $(ODIR), $(SEARCHAF)))
$(BDIR)search-adapter-file: $(searchafobj)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

search-adapter-file: $(BDIR)search-adapter-file


EXEC=filter_reads_with_n search-mgi-adapters trim-mgi-adapters test search-adapter-file
all: $(addprefix $(BDIR), $(EXEC))

.PHONY: clean

clean:
	rm -fr bin/ obj/


