FROM alpine AS build

RUN apk update && apk add g++ gcc musl-dev cmake hdf5-dev make git gfortran python3 && ls /usr/lib/

COPY . /achilles_src

RUN cmake -S achilles_src -B achilles && cmake --build ../achilles -j

FROM alpine AS main

RUN apk update && apk add --no-cache hdf5-dev
COPY --from=build /achilles /achilles
WORKDIR /achilles
ENTRYPOINT ["./bin/achilles"]
CMD ["run.yml"]
