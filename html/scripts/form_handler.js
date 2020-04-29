function openAdvanced() {
    var x = document.getElementById("advanced");
    if (x.style.display === "none") {
        x.style.display = "block";
    } else {
        x.style.display = "none";
    }
}


function validate_read_files() {
    var data = document.getElementById("data").value;

    if (!(data.endsWith('.zip') || data.endsWith('.gz'))){
        alert("Chosen file is illegal. Only zip/tar.gz formats are allowed.");
        return false;
    }
    return true;
}

// function onSubmit(response) {
//     if (validate_read_files() === true) {
//         document.getElementById("microbializer_form").submit();
//     }
// }