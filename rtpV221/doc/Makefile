
all: rtpspec.pdf

rtpspec.pdf: rtpspec.tex
	pdflatex rtpspec.tex -interaction nonstopmode || \
	rm rtpspec.pdf 2> /dev/null || true

show: rtpspec.pdf
	evince rtpspec.pdf > /dev/null 2>&1 &

clean:
	rm *.log *.aux 2> /dev/null || true

