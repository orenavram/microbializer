function openAdvanced() {
    var x = document.getElementById("advanced");
    if (x.style.display === "none") {
        x.style.display = "block";
    } else {
        x.style.display = "none";
    }
}

function validate_lib_file(){
    var alternative_lib = document.getElementById("alternative_lib").value;
    if (alternative_lib != '' && !(alternative_lib.toLowerCase().endsWith('json'))) {
        alternative_lib_file_name = alternative_lib.replace(/^.*[\\\/]/, '');
        alert("Reference library "+alternative_lib_file_name+" is illegal. Only json format is allowed.");
        return false;
    }
    return true;
}
function validate_read_files() {
    var R1s = [document.getElementsByName("run1_R1")[0].value,
        document.getElementsByName("run2_R1")[0].value,
        document.getElementsByName("run3_R1")[0].value,
        document.getElementsByName("run4_R1")[0].value,
        document.getElementsByName("run5_R1")[0].value,
        document.getElementsByName("run6_R1")[0].value];
    var R2s = [document.getElementsByName("run1_R2")[0].value,
        document.getElementsByName("run2_R2")[0].value,
        document.getElementsByName("run3_R2")[0].value,
        document.getElementsByName("run4_R2")[0].value,
        document.getElementsByName("run5_R2")[0].value,
        document.getElementsByName("run6_R2")[0].value];

    for (i = 0; i < R1s.length; i++) {
        if (R1s[i] == '' && R2s[i] == ''){

        } else if (R1s[i] == '' && R2s[i] != '') {
            alert("R1 file of Replicate "+String(i+1)+" is missing (either you upload both of R1 and R2 or none).");
            return false;
        } else if (R1s[i] != '' && R2s[i] == '') {
            alert("R2 file of Replicate "+String(i+1)+" is missing (either you upload both of R1 and R2 or none).");
            return false;
        } else {
            if (!(R1s[i].endsWith('fastq') || R1s[i].endsWith('gz'))) {
                alert("R1 file Replicate "+String(i+1)+" is illegal. Only fastq and gzip formats are allowed.");
                return false;
            }
            if (!(R2s[i].endsWith('fastq') || R2s[i].endsWith('gz'))) {
                alert("R2 file Replicate "+String(i+1)+" is illegal. Only fastq and gzip formats are allowed.");
                return false;
            }
        }
    }
    return true;
}

function setDefaultMassSpecSeq(){
    var MMU = document.getElementsByName("MMU")[0].value;
    if (MMU === "mouse"){
        document.getElementById("mass_spec_seq").value = "AK";
    } else {
        document.getElementById("mass_spec_seq").value = "ASTK";
    }
}