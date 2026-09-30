FROM ubuntu:24.04 AS build
RUN apt-get update && apt-get install -y --no-install-recommends g++ make python3 python3-pip python3-venv \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY src ./src
COPY Makefile ./
RUN make

FROM ubuntu:24.04
RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-flask \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY --from=build /app/bin/bharatsolve ./bin/bharatsolve
COPY ui ./ui
COPY data ./data
COPY demo ./demo
ENV PORT=5055
EXPOSE 5055
CMD ["python3", "ui/server.py"]
