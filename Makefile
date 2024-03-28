.PHONY: help readme tar_misc

help:
	@echo "make help"
	@echo "make readme"
	@echo "make tar_misc"

readme:
	@gh markdown-preview

tar_misc:
	tar czf misc.tar.gz --exclude="*.DS_Store*" misc
