# BioPython & Next.JS project

```sh
docker run --rm -d --name bioseq2 -p 80:80 registry.gitlab.com/mamorosdev/m14-uf2-bioseq-2
```

FAQ's. How to install docker in Linux:

```sh
sudo apt update
sudo apt -y install docker-compose
sudo usermod -aG docker ${USER}
```

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

# Navigation component.

https://medium.com/@a.pirus/how-to-create-a-responsive-navigation-bar-in-next-js-13-e5540789a017
