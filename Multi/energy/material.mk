$(folder)/$(mode)_%.root : $(folder)/cevns_%.xml
	@mkdir -p $(dir $@)
	(export ISOTOPEFILE=$(relicssim)/macro/isotopes/$(ISOTOPE).mac; \
	 cd $(relicssim)/build && \
	 ./BambooMC -c $< -m $(relicssim)/macro/components/$(COMPONENT).mac -n $(eventN) -o $@) \
	 > $@.log 2>&1
