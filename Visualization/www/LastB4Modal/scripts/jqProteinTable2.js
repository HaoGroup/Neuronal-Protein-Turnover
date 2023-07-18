function rdVal(value,digits) { return (value>0 || value<0) ? Math.round(value*(10**digits))/(10**digits) : ''; }
function dplyFixedDigits(value,digits) { 
  // assume value is numeric
  let txt=value.toString();
  let arr=txt.split(".");
  if (arr.length==1) { txt+='.'; arr1len=0; } else { arr1len=arr[1].length; } // some are integers by itself, without decimal
  return txt += Array(digits-arr1len+1).join('0');
}
function roundVal(value,digits) { let val = rdVal(value,digits); return dplyFixedDigits(val,digits); }

function applyFancyTable() {
  proteinTable = $("table#tbProtein").fancyTable({
    sortColumn:1, // columns are 0, 1, 2, ... here
    pagination: true,
    perPage: 50,
    globalSearch: true,
    globalSearchExcludeColumns: [4,5] // columns are 1, 2, 3, ... here. exclude t1/2 and chart columns
  });
  return null;
}

function proteinTableCreate() {
  var tableCreated = $.Deferred();

  // load data
  $.getJSON('data/ProteinRecords.json', function(json) {  
    proteinData = json;
    // colNames = Object.keys(proteinData[0]); // ['prtn', 'gene', 'chart', 'desc', 'peptides'] // 2023-06-12
    // const colNames = ['prtn', 'gene', 'desc','proteinT12', 'chart', 'Npeptides', 'peptides'];
    const colNames = ['prtn', 'gene', 'desc','proteinT12', 'chart', 'peptides'];
    // const colHeaders = {'prtn': 'Protein',  'gene': 'Gene', 'desc':'Description', 'proteinT12':'t<sub>&frac12;</sub> (d)', 'chart': 'Chart', 'Npeptides':'&nbsp;&nbsp;&nbsp;&nbsp;# of<br/>good / all<br/>peptides', 'peptides':'Peptides' } ;
    const colHeaders = {'prtn': 'Protein',  'gene': 'Gene', 'desc':'Description', 'proteinT12':'t<sub>&frac12;</sub> (d)', 'chart': 'Chart', 'peptides':'good / all<br/>peptides' } ;

    // Create Protein main table
    var table = document.getElementById("tbProtein");

    // build table header
    let tbhead = document.createElement("thead");
    let rowhead = document.createElement("tr");
    colNames.forEach( n => { 
      let th = document.createElement('th'); 
      th.className = n; 
      th.setAttribute( 'title', "Click to sort table by this column" );
      if ( n=='proteinT12' || n=='peptides' ) th.setAttribute("data-sortas", "numeric"); // numeric columns need this
      // if ( n=='proteinT12' ) th.setAttribute("data-sortas", "numeric"); // numeric columns need this
      if ( n=="chart" ) th.style.textAlign = "center"; // the check marks are aligned right somehow.
      th.innerHTML = colHeaders[n]; 
      rowhead.appendChild(th); 
    });
    tbhead.appendChild(rowhead);
    table.appendChild(tbhead);

    // build table rows
    let tbbody = document.createElement("tbody");
    proteinData.forEach( proteinrow => {
      const gene = proteinrow['gene']?proteinrow['gene']:proteinrow['prtn'];
      let row = document.createElement("tr");
      row.id = gene;
      tbbody.appendChild(row);
      const chartsN = proteinrow['chart'];
      const t12value = proteinrow['proteinT12'];
      const t12estvalue = proteinrow['proteinT12est'];
      const peptideslen = proteinrow["peptides"]["Peptide"].length;

      // add cells
      colNames.forEach( colname => {
        let cell = document.createElement('td');
        cell.className = colname;
        if (colname == "prtn" || colname == "gene") { 
          cell.innerHTML = proteinrow[colname];
          if (chartsN) {
            if (colname ==  "gene") { cell.style.color = "blue"; }
            cell.setAttribute( 'title', "Click to load charts" );
            cell.setAttribute( 'onclick', 'showProteinChart(this)');
          }
        // } else if (colname == "Npeptides") { 
        //   cell.innerHTML = chartsN+" / "+peptideslen; 
        //   cell.style.color = "blue";
        //   thistitle = "Click to toggle show/hide of peptides.";
        //   cell.setAttribute( 'title', thistitle );
        //   cell.setAttribute( 'onclick', 'pepListToggle(this)');
        } else if (colname == "peptides") {
          // first div to display N-good/N-all
          let pepNdiv = document.createElement("div");
          pepNdiv.style.display = "block"; // show all initially
          pepNdiv.setAttribute("class", "peptideN");
          pepNdiv.style.color = "blue";
          pepNdiv.innerHTML = chartsN+" / "+peptideslen; 
          cell.appendChild(pepNdiv);
          // second div to display list of peptides
          let peplistdiv = document.createElement('div');

          // peplistdiv.style.fontSize="small";
          peplistdiv.style.display = "none"; // hide all initially
          const thisrow = proteinrow[colname];

          // let ul = document.createElement('ul');
          // for (let i=0; i<thisrow['Peptide'].length; i++) {
          //   let li = document.createElement("li");
          //   thissupporti = thisrow['support'][i];
          //   thisr2i = thisrow['r2_CFit'][i];
          //   thist12i = thisrow['t12_CFit'][i];
          //   goodpep = (thissupporti>2 && thisr2i>0.8)?true:false;
          //   li.innerHTML = thisrow['Peptide'][i] + ' ('+thissupporti+') t<sub>&frac12;</sub>: ' + roundVal(thist12i,2)+ 'd, r<sup>2</sup>: '+ roundVal(thisr2i,3) ;
          //   if (goodpep) { li.style.fontWeight="bold"; } else { li.style.fontStyle="italic"; } 
          //   ul.appendChild(li);
          // }
          // peplistdiv.appendChild(ul);

          // use a table to display list of peptides
          let tb2 = document.createElement("table");
          let tb2head = document.createElement("thead");
          const tb2ColHeaders = { 'peptide':'Peptide', 'support':'Support', 't12':'t<sub>&frac12;</sub>', 'r2':'r-sq' }
          let rowhead = document.createElement("tr");
          Object.keys(tb2ColHeaders).forEach( n => { 
            let th2 = document.createElement('th'); 
            th2.className = 'peptb_'+n; 
            // th.setAttribute( 'title', "Click to sort table by this column" );
            // if ( n=="chart" ) th.style.textAlign = "center"; // the check marks are aligned right somehow.
            th2.innerHTML = tb2ColHeaders[n]; 
            rowhead.appendChild(th2); 
          });
          tb2head.appendChild(rowhead);
          tb2.appendChild(tb2head);
          let tb2body = document.createElement("tbody");
          // let tr = document.createElement("tr");
          for (let i=0; i<thisrow['Peptide'].length; i++) {
            let tr = document.createElement("tr");
            let tdpep = document.createElement("td");
            let tdsupport = document.createElement("td");
            let tdt12 = document.createElement("td");
            let tdr2 = document.createElement("td");
            let thissupporti = thisrow['support'][i];
            let thisr2i = thisrow['r2_CFit'][i];
            let thist12i = thisrow['t12_CFit'][i];

            tdpep.innerHTML = thisrow["Peptide"][i]; 
            tdsupport.innerHTML = thissupporti;
            tdt12.innerHTML = thist12i;
            tdr2.innerHTML = thisr2i;

            goodpep = (thissupporti>2 && thisr2i>0.8)?true:false;
            if (goodpep) { tr.style.fontWeight="bold"; } else { tr.style.fontStyle="italic"; } 
            tr.appendChild(tdpep);
            tr.appendChild(tdsupport);
            tr.appendChild(tdt12);
            tr.appendChild(tdr2);
            tb2body.appendChild(tr);
          }
          tb2.appendChild(tb2body);
          // debugger;
          peplistdiv.appendChild(tb2);

          peplistdiv.setAttribute("class", "peptideList");
          cell.appendChild(peplistdiv);
          //
          thistitle = "Click to toggle show/hide of peptides. Bold face indicates the 'good' peptides.";
          cell.setAttribute( 'title', thistitle );
          cell.setAttribute("class", "peptideCell");
          cell.setAttribute( 'onclick', 'pepListToggle(this)');
          // cell.set
        } else if (colname == "chart") {
          // cell.innerHTML = (proteinrow[colname]>0) ? '&#x2713;' : '&#x274C;'; // 2713/2714
          if (chartsN) {
            cell.style.color = "green";
            cell.innerHTML = '&#x2713;' ;  // '&#x2713;', '&#x2714;'  
          } else {
            cell.style.color = "red";
            cell.innerHTML = '&#x2717;' ; //  '&#x2717;', '&#x274C;
          }
        } else if (colname == 'proteinT12') {
          if (t12value) {
            cell.innerHTML = roundVal(t12value,2);//
            cell.style.fontWeight = "bold";
            thiscelltitle = "Average (harmonic mean) of the 'good' peptide half lives. If ALL peptides are used, the estimated average (harmonic mean) is "+roundVal(t12estvalue,2)+"d" ;
          } else {
            cell.innerHTML = roundVal(t12estvalue,2);//
            cell.style.fontStyle = "italic";
            thiscelltitle = "No 'good' peptide data series available. Using ALL peptide series, the estimated average half life (harmonic mean) is "+roundVal(t12estvalue,2)+"d" ;
          }
          cell.setAttribute( 'title', thiscelltitle);
        } else { cell.innerHTML = proteinrow[colname]; } // others like 'desc'
        row.appendChild(cell);
      });    
    });
    // dataLoaded.resolve();
    table.appendChild(tbbody);
    let div = $("div#divProteinTable")[0];
    div.innerHTML = null; // remove "Loading table ..."
    div.appendChild(table); // div is an array

    // browser=browserChk(true);
    // if ( !( browser['isChrome'] || browser['isFirefox'] || browser['isIE'] || browser['isEdge'] || browser['isBlink'] )) {
    //   console.log('checked: not chrome, firefox, IE, edge, nor Blink');
    // };

    tableCreated.resolve();
  });

  return tableCreated.promise();
}

function browserChk(log=false) {
  // Safari, for example, complains about Maximum call stack size exceeded with FancyTable iterating thru 10k of data rows. 
  // check broswer: https://stackoverflow.com/questions/9847580/how-to-detect-safari-chrome-ie-firefox-and-opera-browsers
  var isOpera = (!!window.opr && !!opr.addons) || !!window.opera || navigator.userAgent.indexOf(' OPR/') >= 0; // Opera 8.0+
  var isFirefox = typeof InstallTrigger !== 'undefined'; // Firefox 1.0+
  var isSafari = /constructor/i.test(window.HTMLElement) || (function (p) { return p.toString() === "[object SafariRemoteNotification]"; })(!window['safari'] || (typeof safari !== 'undefined' && window['safari'].pushNotification)); // Safari 3.0+ "[object HTMLElementConstructor]" 
  var isIE = /*@cc_on!@*/false || !!document.documentMode; // Internet Explorer 6-11
  var isEdge = !isIE && !!window.StyleMedia; // Edge 20+
  var isChrome = !!window.chrome; // !!window.chrome && (!!window.chrome.webstore || !!window.chrome.runtime); // Chrome 1 - 79
  var isEdgeChromium = isChrome && (navigator.userAgent.indexOf("Edg") != -1); // Edge (based on chromium) detection
  var isBlink = (isChrome || isOpera) && !!window.CSS;  // Blink engine detection

  if (log) { console.log(' isOpera: '+isOpera+'\n isFirefox: '+isFirefox+'\n isSafari: '+isSafari+'\n isIE: '+isIE+'\n isEdge: '+isEdge+'\n isChrome: '+isChrome+'\n isEdgeChromium: '+isEdgeChromium+'\n isBlink: '+isBlink); }

  return {isOpera:isOpera, isFirefox:isFirefox, isSafari:isSafari, isIE:isIE, isEdge:isEdge, isChrome:isChrome, isEdgeChromium:isEdgeChromium, isBlink:isBlink}
}

function showProteinChart(elm) {
  var row = elm.parentNode;
  const gene = row.id; // event.target.attributes['proteinlink'].value;

  // set focus 
  const position =  $("div#proteinChart iframe").position();
  scroll(0,position.top);

  var proteinframe = $("div#proteinChart iframe")[0];
  var peptideframe = $("div#peptideChart iframe")[0];
  const proteinPage = "/charts/proteinLevel/RelAbundance_Gene-"+gene+".html";
  const peptidePage = "/charts/peptideLevel/RelAbundance_Gene-"+gene+"-peptides.html";

  proteinframe.setAttribute('src', proteinPage);
  peptideframe.setAttribute('src', peptidePage);

  // change the pop-up links
  var proteinpopup = $("div#proteinChart a:first")[0];
  var peptidepopup = $("div#peptideChart a:first")[0];
  proteinpopup.setAttribute('href', proteinPage);
  peptidepopup.setAttribute('href', peptidePage);

  return null;
}

function pepListToggle(elm) {
  // var row = elm.parentNode;
  var row = elm;
  var tdPepList = row.querySelector("div.peptideList");
  var tdPepN = row.querySelector("div.peptideN");
  var curStatus = tdPepList.style.display;
  if (curStatus==="none") { 
    $("td > div.peptideN").show();
    $("td > div.peptideList").hide();
    tdPepList.style.display = "block"; 
    tdPepN.style.display="none";
  } else { // curStatus=="block"
    $("td > div.peptideN").show();
    $("td > div.peptideList").hide();
  }
  return null;
}

function showFinalTable() {
  browser=browserChk(false);
  if ( !( browser['isChrome'] || browser['isFirefox'] || browser['isIE'] || browser['isEdge'] || browser['isBlink'] )) {
    console.log('checked: not chrome, firefox, IE, edge, nor Blink');
    return null;
  };
  
  let t1 = document.getElementById('loadingMsg');
  t1.innerHTML = "&nbsp;&nbsp;Sortable columns (Click on a column header to sort)";
  let t2 = document.getElementById('divProteinTable');
  t2.style.display = 'block';
  return null;
}
