process printit{

    input:
    val x
    
    script:
    """
    echo ${x}
    """
}