SHELL:=bash
fileN:=$(FILEN)
eventN:=$(EVENTN)
folder:=$(FOLDER)
geo:=$(GEO)
suffix:=$(SUFFIX)
mode:=$(PMODE)
relicssim:=$(RELICSSIM)
analysis:=$(ANALYSIS)

format:=07g
filenum:=$(shell seq -f '%$(format)$(suffix)' 1 $(fileN))
json:=$(filenum:%=$(folder)/params_%.json)
combine:=$(folder)/$(mode)$(suffix).root
collection:=$(folder)/events$(suffix).h5

.PHONY : all

all : $(collection)

$(folder)/params_%.json :
	@mkdir -p $(dir $@)
	cat $(relicssim)/config/geo_params.json $(geo) > $@

$(folder)/cevns_%.xml : $(folder)/params_%.json
	@mkdir -p $(dir $@)
	python3 $(relicssim)/scripts/relics_xml.py $< -o $@ --run_number $(shell echo $* | cut -d'_' -f 1) --gen $(mode)

include $(mode).mk

$(folder)/events_%.h5 : $(folder)/$(mode)_%.root
	@mkdir -p $(dir $@)
	OMP_NUM_THREADS=2 python3 $(analysis)/Sorting.py --InputFile $< --OutputFile $@ > $@.log 2>&1

$(collection) : $(patsubst $(folder)/params_%.json,$(folder)/events_%.h5,$(json))
	@mkdir -p $(dir $@)
	OMP_NUM_THREADS=2 python3 $(analysis)/Combine.py --InputFiles $(folder)/events_%$(format)$(suffix).h5 --OutputFile $@ -N $(fileN) > $@.log 2>&1

$(combine) : $(patsubst $(folder)/params_%.json,$(folder)/$(mode)_%.root,$(json))
	@mkdir -p $(dir $@)
	hadd $(combine) $^

clean :
	rm -rf $(folder)/*$(suffix)*

.DELETE_ON_ERROR: 

.SECONDARY:
