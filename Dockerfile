FROM alpine as build

RUN apk update && apk add g++ gcc musl-dev cmake hdf5-dev make git gfortran && ls /usr/lib/

COPY . /achilles_src

RUN cmake -S achilles_src -B achilles -DENABLE_BSM=OFF && cmake --build ../achilles -j

FROM alpine as main

RUN apk update && apk add --no-cache hdf5-dev
COPY --from=build /achilles /achilles
WORKDIR /achilles
ENTRYPOINT ["./bin/achilles"]
CMD ["run.yml"]
