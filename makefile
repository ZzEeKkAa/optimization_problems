HOST := 'm3'
remoute-build:
	GOOS=linux GOARCH=amd64 go build -o ./lab2/lab2 ./lab2/
	scp ./lab2/lab2 $(HOST):~/lab2
	rm ./lab2/lab2
remoute-run:
	ssh $(HOST) './lab2 -threads=32 -iterations=10000'
	scp $(HOST):img.gif ./img.gif