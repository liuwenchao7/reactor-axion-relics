$(folder)/$(mode)_%.root : $(folder)/cevns_%.xml
	@mkdir -p $(dir $@)
	$(relicssim)/build/BambooMC -c $< -n $(eventN) -o $@ > $@.log 2>&1
