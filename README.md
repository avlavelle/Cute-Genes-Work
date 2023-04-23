# Gene Party 3.0: Parsing HGNC Data and Using Kubernetes

This project is compiled of one script that builds a Flask application and provides different routes for the user to query that return different pieces of information about the HGNC data, using the HGNC data from HUGO.The project also includes a Dockerfile that allows any user to pull the image from Docker Hub and run an instance of the image into a container on their personal machine, as well as a Docker-compose file that allows a user to start the container easily and succinctly.

The Flask application routes are important because the HGNC data contains loads of information that may not be interesting to the user. The Flask application and the routes allow users to quickly query the information that they may be looking for without having to comb through the data manually. The application also utilizes the Redis database, making it possible to store and come back to the data. The Dockerfile allows any user to have access to these capabilities through running an instance of the image on their own machine, plus the Docker compose file automates the deployment of the application, so the container is easy to start. 

The project also includes five .yml files that include a Redis service, deployment, and pvc, as well as a Flask service and deployment. These files allow the project to be deployed into Kubernetes.

## Accessing and Describing the Data
 
The HGNC data is loaded in from this link, https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/json/hgnc\_complete\_set.json , where the data is in json format. The application uses requests to access the data, and the POST route puts the data into the Redis database. The data contains lots of unique and interesting information about every hgnc\_id that the HGNC knows exists.

## Flask App and Its Routes

The Flask application contains a function that loads the data into the Redis database and different routes that return all of the information, only hgnc\_id's, or all of the information about a specific hgnc\_id. There is also a route that clears the data from the database.

## Pulling the Image and Running the App on a VM

Install the project by cloning the repository. 

Use ```docker pull avlavelle/gene_api``` to pull a copy of the container.
Then, ```docker-compose up --build``` will get the container running using the compose file, build the image, and bind the appropriate port inside the container to the appropriate port on the host.

In a separate window, you can use ``` curl localhost:5000/<route> ``` to call the routes.

## Building a New Image from the Dockerfile

In order to build a new image from the Dockerfile, use the same ```docker pull``` command from above. 
Then, use ```docker-compose build``` to build the image and use ```docker images``` to check that a copy of the image has been built.

## Using and Deploying into Kubernetes

Clone the repository and pull the image in a Kubernetes cluster using instructions from above. No changes should be required for the Kubernetes files or the Python script because an environment variable should dynamically set the Redis host IP. 

Users must then ```kubectl apply -f <file_name.yml>``` for each service, deployment, and pvc. 

Then, users should exec into the debug pod. Using ```kubectl get pods```, the debug bod name should become available. Exec into this pod using ```kubectl exec -it <debug_pod_name> --/bin/bash```.

Users should then be able to ```curl avlav-test-flask-service:5000/<route>``` all of the routes from within the debug pod.

## ```/image``` Route

The ```/image -X POST``` route reads the date each unique gene entry was first approved from the database, tabulates how many genes were approved each year, and creates a bar graph of the data, which it writes into another database.

The ```/image``` route is the default GET request. It returns the plot to the user, so when the route is called, the user must name their plot in the ```curl``` request. For example, ```curl localhost:5000/image>>myimage.png```. After that, the image should be created and available within the repository on the VM. Two ```scp``` actions are then required to see the image.

The first takes place on the VM:
```
scp ./ <image_name.png> username@address.edu:./
```
The commandline should return confirmation of the ```scp```. 

On the user's local machine in the intended folder or repository:
```
scp username@address.edu:./<image_name.png> ./
```
The commandline should return confirmation of the ```scp```, and the user should be able to access the plot image from their file explorer now.

The ```/image -X DELETE``` route deletes the image from the database.

## Example Queries and Interpretation of Results

For the ```/data``` POST route which puts the data into Redis:
```
Data loaded
```

For the ```/data``` DELETE route which deletes the data in Redis:
```
Data deleted, there are 0 keys in the db
```

For the ```/data``` GET route which returns all the data from Redis:
```
[{"version": ..., "alias_symbol": [...], "data_approved_reserved": ...}...]
```

For the ```/genes``` route which returns a json-formatted list of all hgnc_ids:
```
[..., "HGNC:13195", ..., "HGNC:24523"]
```

For the ```/genes/<hgnc_id>``` route which returns all the data associated with <hgnc_id>:
```
{["hgnc_id": "HGNC:24523", "location": ..., ...]}
```

For the ```/image``` POST route which creates a plot image and loads it into Redis:
```
Image created
```

For the ```/image``` GET route which retrieves the plot:
The user should follow the instructions above to retreive the image.
The return of this route will not be visible to the user. 

For the ```/image``` DELETE route which deletes the plot from the database:
```
Image deleted, there are 0 keys in the db
```
