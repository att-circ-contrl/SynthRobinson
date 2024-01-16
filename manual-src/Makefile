# Data synthesis using the Robinson/Freyer/Hindriks model - Documentation.
# Written by Christopher Thomas.

REFNAME=rfh-libref
REFBASE=ref

GUIDENAME=rfh-guide
GUIDEBASE=guide


GUIDEBITS=\
	$(GUIDENAME).tex		\
	$(GUIDEBASE)-over.tex		\
	$(GUIDEBASE)-model.tex		\
	$(GUIDEBASE)-howto.tex		\
	$(GUIDEBASE)-math.tex

REFBITS=\
	$(REFNAME).tex

# No overview for now; the not much to say in one.
#	$(REFBASE)-over.tex

REFAUTO=\
	$(REFBASE)-funcs.tex		\
	$(REFBASE)-examples.tex		\
	$(REFBASE)-notes.tex


default: force

guide: $(GUIDENAME).pdf
libref: $(REFNAME).pdf

force: clean guide libref install
	make autoclean

clean: tclean autoclean
	rm -f $(GUIDENAME).pdf
	rm -f $(REFNAME).pdf

autoclean:
	rm -f $(REFAUTO)

tclean:
	rm -f $(GUIDENAME)*.log $(GUIDENAME)*.aux $(GUIDENAME)*.toc
	rm -f $(REFNAME)*.log $(REFNAME)*.aux $(REFNAME)*.toc
	rm -f texput.*

install:
	mv *pdf ../manual


$(GUIDENAME).pdf: $(GUIDEBITS)
	pdflatex $(GUIDENAME).tex
	pdflatex $(GUIDENAME).tex
	make tclean

refauto-bits: refauto-funcs refauto-examples refauto-notes

LIBDIR=../library
EXAMPLEDIR=../examples

refauto-funcs:
	./copy-code.pl --outfile=$(REFBASE)-funcs.tex --matlabdoc \
		"--chapter=Function Reference" --label=sect-funcs \
		$(LIBDIR)/*.m

refauto-examples:
	./copy-code.pl --outfile=$(REFBASE)-examples.tex --verbatim \
		"--chapter=Sample Code" --label=sect-examples \
		$(EXAMPLEDIR)/*/do*.m

refauto-notes:
	./copy-code.pl --outfile=$(REFBASE)-notes.tex --verbatim \
		"--chapter=Data Structures and Additional Notes" \
		--label=sect-notes `find $(LIBDIR)|grep txt`

$(REFNAME).pdf: $(REFBITS) refauto-bits
	pdflatex $(REFNAME).tex
	pdflatex $(REFNAME).tex
	make tclean

#
# This is the end of the file.
