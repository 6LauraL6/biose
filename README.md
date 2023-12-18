# BioPython & Next.JS project

```sh
docker run --rm -d --name kekule -p 80:80 registry.gitlab.com/xtec/bio-sequence
```

Imatge també disponible a:
```sh
docker run --rm -d --name biosequence2 -p 80:80 registry.gitlab.com/mamorosdev/m14-uf2-bioseq-2
```

## Deployed at:

https://dawbio-m14-flask-pt2.azurewebsites.net/


## Develop

Install dependencies:

```sh
./init.sh
source .bashrc
```


## Run:

```sh
npm run dev
```

or:

```sh
npm run next-dev
npm run flask-dev
```