-- This has to be run in a psql shell. Connect with
-- psql -h {hostname} -U {user} -p {port} {db name}
-- Then copy and paste line by line to get the data into the database. 

\COPY soil.model_error FROM '~/git/soil-retention-curves/data/for_db/model_error.csv' DELIMITER ',' CSV HEADER
\COPY soil.parameters FROM '~/git/soil-retention-curves/data/for_db/parameters.csv' DELIMITER ',' CSV HEADER
\COPY soil.raw_soil_data FROM '~/git/soil-retention-curves/data/for_db/vwc_kpa_observations.csv' DELIMITER ',' CSV HEADER
