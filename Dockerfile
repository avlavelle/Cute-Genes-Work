FROM python:3.8.10

RUN pip install Flask==2.2.0 redis requests
RUN pip install requests==2.22.0
RUN pip install collection==0.1.6
RUN pip install numpy==1.24.2
RUN pip install matplotlib==3.7.1


ADD gene_api.py /gene_api.py

CMD ["python", "gene_api.py"]
