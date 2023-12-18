# Bio Kekule

Provem de fer un canvi i a veure si funciona CI

```sh
docker run --rm -d --name kekule -p 80:80 registry.gitlab.com/xtec/bio-sequence
```

## Develop

Install dependencies:

```sh
./init.sh
source .bashrc
```


Run:

```sh
npm run dev
```

or:

```sh
npm run next-dev
npm run flask-dev
```