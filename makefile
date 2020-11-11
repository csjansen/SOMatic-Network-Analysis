FILES = MotifZscore MakeNetwork

build:
	for dir in $(FILES); do \
		$(MAKE) -C $$dir; \
	done
clean:
	for dir in $(FILES); do \
		$(MAKE) -C $$dir clean; \
	done
