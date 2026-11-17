.PHONY: pdf check integrity clean

LATEX = pdflatex -interaction=nonstopmode -halt-on-error

pdf:
	cd sixth-moment && $(LATEX) moments.tex
	cd sixth-moment && $(LATEX) moments.tex
	cd paper && $(LATEX) moments.tex
	cd paper && $(LATEX) moments.tex
	cd supplement && $(LATEX) moments.tex
	cd supplement && $(LATEX) moments.tex
	cd summary && $(LATEX) summary.tex
	cd summary && $(LATEX) summary.tex

integrity:
	python3 checks/check_manifest.py

check: integrity
	python3 checks/check_ten_terms.py
	python3 checks/derive_ten_terms.py
	python3 checks/check_finite_derivative.py
	python3 checks/exponents.py
	python3 checks/identity_checks.py
	python3 checks/elementary_pair_checks.py

clean:
	rm -f paper/*.aux paper/*.log paper/*.out sixth-moment/*.aux sixth-moment/*.log sixth-moment/*.out supplement/*.aux supplement/*.log supplement/*.out summary/*.aux summary/*.log summary/*.out
