# Build stage
FROM rust:alpine AS builder

RUN apk add --no-cache musl-dev pkgconfig openssl-dev openssl-libs-static cmake make g++ perl linux-headers

WORKDIR /app
COPY Cargo.toml Cargo.lock ./
COPY src ./src

RUN cargo build --release

# Runtime stage - bare metal
FROM scratch

COPY --from=builder /app/target/release/fasterp /fasterp

ENTRYPOINT ["/fasterp"]
