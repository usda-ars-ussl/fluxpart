TUTORBUILD    = tutor.py
RESULTSDIR    = results
TXTDIR        = ../text
FIGDIR        = ../images

.DEFAULT_GOAL := tutor

.PHONY: clean
clean:
	rm -rf $(RESULTSDIR)/* $(TXTDIR)/* $(FIGDIR)/*

.PHONY: tutor
tutor:
	python $(TUTORBUILD)
	@echo
	@echo "Tutorial build finished."
