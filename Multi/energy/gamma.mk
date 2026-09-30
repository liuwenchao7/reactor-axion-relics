$(folder)/$(mode)_%.root : $(folder)/cevns_%.xml
	@mkdir -p $(dir $@)
	GAMMAHIST=$(relicssim)/data/gamma$(REACTOR).txt $(relicssim)/build/BambooMC -c $< -m $(relicssim)/macro/gamma$(TOPSIDE).mac -n $(eventN) -o $@ > $@.log 2>&1
