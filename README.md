# metroman-docker

MetroMan docker container hosted on AWS ECR used to create containers in the cloud.

***NOTE: Confluence version input operations currently under development***

# publishing disclaimer

Please **do not share this** repository, algorithm code contained within, or containers with others. 

***Also note that all Docker containers are currently hosted in AWS. You may build the containers using the Dockerfiles located in each algorithm's directory.***

# Container

## Arguments

Each container takes the name of a reach file as an argument. This file contains a list of reaches or sets of reaches where each reach identifier or set is on one line.

### Prediagnostics
- EFS Data is mounted to the container at `/mnt/data`
- (Note data is overwritten)

## AWS Environment variable

AWS_BATCH_JOB_ARRAY_INDEX is used to determine the reach identifier to process for all containers. The index value is used to retrieve a line in a file that lists all reach identifiers to be processed.