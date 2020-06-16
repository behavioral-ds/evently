# generate a twitter-snowflake id, based on
# https://github.com/twitter/snowflake/blob/master/src/main/scala/com/twitter/service/snowflake/IdWorker.scala
make_snowflake <- function(timestamp_ms, datacenter_id, worker_id, sequence_id, twepoch = as.integer64('1288834974657')) {
  check_required_packages('bit64')
  base <- as.integer64(2)
  datacenter_id_bits <- 5
  worker_id_bits <- 5
  sequence_id_bits <- 12
  max_datacenter_id <- 1 * base^datacenter_id_bits
  max_worker_id <- 1 * base^worker_id_bits
  max_sequence_id <- 1 * base^sequence_id_bits
  max_timestamp <- base^(64 - datacenter_id_bits - worker_id_bits - sequence_id_bits)

  stopifnot(is.character(timestamp_ms) || is.integer64(timestamp_ms))
  timestamp_ms <- as.integer64(timestamp_ms)
  sid <- ((timestamp_ms - twepoch) %% max_timestamp) * base^datacenter_id_bits * base^worker_id_bits * base^sequence_id_bits
  sid <- sid + (datacenter_id %% max_datacenter_id) * base^worker_id_bits * base^sequence_id_bits
  sid <- sid + (worker_id %% max_worker_id) * base^sequence_id_bits
  sid <- sid + sequence_id %% max_sequence_id
  sid
}

# inversely transform a snowflake id back to its components.
melt_snowflake <- function(snowflake_id, twepoch = as.integer64('1288834974657')) {
  check_required_packages('bit64')
  base <- as.integer64(2)
  datacenter_id_bits <- 5
  worker_id_bits <- 5
  sequence_id_bits <- 12
  max_datacenter_id <- 1 * base^datacenter_id_bits
  max_worker_id <- 1 * base^worker_id_bits
  max_sequence_id <- 1 * base^sequence_id_bits
  max_timestamp <- base^(64 - datacenter_id_bits - worker_id_bits - sequence_id_bits)

  stopifnot(is.character(snowflake_id) || is.integer64(snowflake_id))
  snowflake_id <- as.integer64(snowflake_id)
  sequence_id <- snowflake_id %% max_sequence_id
  worker_id <- (snowflake_id %/% base^sequence_id_bits) %% max_worker_id
  datacenter_id <- (snowflake_id %/% base^sequence_id_bits %/% base^worker_id_bits) %% max_datacenter_id
  timestamp_ms <- snowflake_id %/% base^sequence_id_bits %/% base^worker_id_bits %/% base^datacenter_id_bits
  timestamp_ms <- timestamp_ms + twepoch

  list(timestamp_ms = timestamp_ms,
       datacenter_id = datacenter_id,
       worker_id = worker_id,
       sequence_id = sequence_id)
}
