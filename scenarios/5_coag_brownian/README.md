How to run the scenario using Docker?

```bash
$ mkdir out
$ docker run -it --rm -v $PWD:/run compdyn/partmc bash -c 'cd /run; /build/partmc example.spec'
```

