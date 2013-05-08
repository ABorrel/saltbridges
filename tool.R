#!/usr/bin/env Rscript

defColor = function (l_name){
	out = c()
	for (element in l_name){
		print (element)
		if (element == "OxAcid"){
			out = append (out, "red")
		}
		else if (element == "ODonAcc"){
			out = append (out, "orange")
		}
		else if (element == "Sulfur"){
			out = append (out, "gold")
		}
		else if (element == "Nhis"){
			out = append (out, "darkblue")
		}
		else if (element == "Nbasic"){
			out = append (out, "blue")
		}
		else if (element == "Carom"){
			out = append (out, "purple")
		}
		else if (element == "OxPep"){
			out = append (out, "coral")
		}
		else if (element == "Ndonnor"){
			out = append (out, "green")
		}
		else if (element == "OxAccept"){
			out = append (out, "lemonchiffon")
		}
		else if (element == "H2O"){
			out = append (out, "cyan")
		}
		else if (element == "CPep"){
			out = append (out, "brown")
		}
		else if (element == "others"){
			out = append (out, "grey")
		}
		else {
			out = append (out, "black")
		}
	}
	names(out) = l_name
	return (out)
}

