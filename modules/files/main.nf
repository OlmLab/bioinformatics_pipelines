process concatenate_files {
    /*
    * This process takes in a list of files and concatenates them into a single file.
    */
    publishDir params.output_dir
    input:
    path file_list
    val output_file

    output:
    path "${output_file}", emit: concatenated_file

    script:
    """
    cat ${file_list} > ${output_file}
    """
}
def tableToDict(file, delimiter = ',') {
    /*
    * This function reads a CSV file and converts it into a dictionary (map) where the keys are the headers
    * and the values are lists of the corresponding column values.
    * @param file: The CSV file to read.
    * @param delimiter: The delimiter used in the CSV file (default is comma).
    * @return: A map where keys are headers and values are lists of column values.
    */



    def result = [:]
    def lines = file.readLines()
    
    // Get all headers
    def headers = lines[0].split(delimiter).collect { it.trim() }
    
    // Initialize empty arrays for each header
    headers.each { header ->
        result[header] = []
    }
    
    // Collect all values for each column
    lines[1..-1].each { line ->
        def values = line.split(delimiter).collect { it.trim() }
        headers.eachWithIndex { header, i ->
            result[header] << values[i]
        }
    }
    
    return result
}
