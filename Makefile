all: docker

FILE=interval_truth_model

## build docker image (requires root rights for docker)
dbuild: Dockerfile
	docker build \
	-t $(FILE) .

## run docker container
docker: dbuild
	echo "open RStudio Server at http://localhost:8787"
	docker run \
	--rm \
	-ti \
	-e DISABLE_AUTH=true \
	-e ROOT=true \
	-e USERID=$(id -u) \
	-e GROUPID=$(id -g) \
	-p 8787:8787 \
	-v $(CURDIR)/sim_results:/home/rstudio/sim_results \
	-v $(CURDIR)/src:/home/rstudio/src \
	$(FILE)
