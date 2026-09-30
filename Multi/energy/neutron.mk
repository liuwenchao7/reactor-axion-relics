$(folder)/$(mode)_%.root : $(folder)/cevns_%.xml
	@mkdir -p $(dir $@)
	NEUTRONHIST=$(relicssim)/data/neutron$(REACTOR).txt $(relicssim)/build/BambooMC -c $< -m $(relicssim)/macro/neutron$(TOPSIDE).mac -n $(eventN) -o $@ > $@.log 2>&1
