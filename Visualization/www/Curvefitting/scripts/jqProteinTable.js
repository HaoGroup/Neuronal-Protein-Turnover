function rdVal(value,digits) { return (value>0 || value<0) ? Math.round(value*(10**digits))/(10**digits) : ""; }
function dplyFixedDigits(value,digits) { 
  // assume value is numeric
  if (!value) return "&nbsp;"; 
  let txt=value.toString();
  let arr=txt.split(".");
  if (arr.length==1) { txt+='.'; arr1len=0; } else { arr1len=arr[1].length; } // some are integers by itself, without decimal
  return txt += Array(digits-arr1len+1).join('0');
}
function roundVal(value,digits) { let val = rdVal(value,digits); return dplyFixedDigits(val,digits); }

function applyFancyTable() {
  proteinTable = $("table#tbProtein").fancyTable({
    sortColumn: 1, // columns are 0, 1, 2, ... here
    pagination: true,
    // perPage: 50,
    globalSearch: true,
    globalSearchExcludeColumns: [4,5] // columns are 1, 2, 3, ... here. exclude t1/2 and chart columns
  });
  return null;
}

function proteinTableCreate() {
  var tableCreated = $.Deferred();

  // load data
  $.getJSON('data/ProteinSummary.json', function(json) {  
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
      const rowjson = JSON.stringify(proteinrow);
      const gene = proteinrow['gene']?proteinrow['gene']:proteinrow['prtn'];
      let row = document.createElement("tr");
      row.id = gene;
      tbbody.appendChild(row);
      // debugger;

      const model = 'CFit'
      const chartsN = proteinrow['chart'];
      const t12value = proteinrow['t12_'+model+'_pass'];
      const t12allvalue = proteinrow['t12_'+model+'_all'];
      const t12bestvalue = proteinrow['t12_'+model+'_best'];
      const peptideslen = proteinrow["peptides"]["Peptide"].length;

      // add cells
      colNames.forEach( colname => {
        let cell = document.createElement('td');
        cell.className = colname;
        if (colname == "prtn" || colname == "gene") { 
          cell.innerHTML = proteinrow[colname];
          if (chartsN) {
            cell.style.cursor = "pointer"; // for both Protein and Gene columns
            if (colname ==  "gene") { cell.style.color = "blue"; }
            // cell.setAttribute( 'title', "Click to load charts" );
            cell.setAttribute( 'onclick', 'showProteinChart(this)');
          }
        } else if (colname == "peptides") {
          let button = document.createElement("button");
          button.setAttribute("class", "btn btn-link");
          button.setAttribute("data-bs-toggle", "modal");
          button.setAttribute("data-bs-target", "#modalPeptideList");
          button.setAttribute("data-bs-peptides", rowjson);
          button.innerHTML = chartsN+" / "+peptideslen; 
          cell.appendChild(button);
          //
          thistitle = "Click to toggle show/hide of peptides. Bold face indicates the 'good' peptides.";
          cell.setAttribute( 'title', thistitle );
          cell.setAttribute("class", "peptideCell");
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
            thiscelltitle = "Average (harmonic mean) of the 'good' peptide half lives. If ALL peptides are used, the estimated average (harmonic mean) is "+roundVal(t12allvalue,2)+"d" ;
          } else {
            cell.innerHTML = roundVal(t12bestvalue,2);//
            cell.style.fontStyle = "italic";
            thiscelltitle = "No 'good' peptide data series available. Using ALL peptide series, the estimated average half life (harmonic mean) is "+roundVal(t12bestvalue,2)+"d" ;
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

    // peptide List modal popup 
    var modalPeptideList = document.getElementById('modalPeptideList')
    modalPeptideList.addEventListener('show.bs.modal', function (event) {
      // var button = event.relatedTarget;  // Button that triggered the modal
      var proteinrow = JSON.parse( event.relatedTarget.getAttribute('data-bs-peptides') ); // Extract info from data-bs-* attributes
      // If necessary, you could initiate an AJAX request here // and then do the updating in a callback.

      const model = 'CFit'
      const prtn = proteinrow['prtn'];
      const gene = proteinrow['gene'];
      const chartsN = proteinrow['chart'];
      const t12value = proteinrow['t12_'+model+'_pass'];
      const t12allvalue = proteinrow['t12_'+model+'_all'];
      const t12bestvalue = proteinrow['t12_'+model+'_best'];
      const peptideslen = proteinrow["peptides"]["Peptide"].length;

      // also refresh charts 
      if (chartsN) showProteinChart(elm=null, genename=gene);

      // set modal title
      var modalTitle = modalPeptideList.querySelector('.modal-title');
      modalTitle.textContent = 'Protein / Gene: ' + proteinrow['prtn'] + ' / ' + proteinrow['gene'];
      // set modal body peptide Good count
      var modalPeptideGdCnt = modalPeptideList.querySelector('.modal-body ul li#modalPeptideGdCnt');
      modalPeptideGdCnt.innerHTML = 'No. of peptides with good data:<span style="background-color:#FEF1C2; font-weight:bold;">&nbsp;&nbsp;'+chartsN+'&nbsp;&nbsp;</span>' ; 
      // modalPeptideGdCnt.setAttribute("class","bg-warning text-dark"); 
      // set modal body peptide All count
      var modalPeptideAllCnt = modalPeptideList.querySelector('.modal-body ul li#modalPeptideAllCnt');
      modalPeptideAllCnt.innerHTML = 'No. of ALL peptides: <i>'+peptideslen+'</i>' ; 
      // set modal body t-half
      var modalProteinT12 = modalPeptideList.querySelector('.modal-body ul li#modalProteinT12');
      modalProteinT12.innerHTML = (t12value) ? 'Average half life (good/all peptides) =<span style="background-color:#FEF1C2; font-weight:bold;">&nbsp;&nbsp;' + roundVal(t12value,2)+' d </span>/ <span style="font-style:italic">' +roundVal(t12allvalue,2) + ' d</span>' : 'Not enough good data from any peptide. Estimated average half life from ALL peptides ≈ <span style="font-style:italic">' + roundVal(t12bestvalue,2) + ' d</span>';
      // set modal body protein description
      var modalProteinDesc = modalPeptideList.querySelector('.modal-body ul li#modalProteinDesc');
      modalProteinDesc.innerHTML = 'Protein Desscription: '+proteinrow['desc'] ; 
      // set peptide list detail table
      var modalPeptideTb = modalPeptideList.querySelector('.modal-body p#modalPeptideTb');

      let pepTable = document.createElement("table"); pepTable.setAttribute("class","table table-hover");

      let pepTablehead = document.createElement("thead");
      const pepTableColHeaders = { 'peptide':'Peptide', 'support':'Support', 't12':'t<sub>&frac12;</sub> (d)', 'r2':'r-squared' }
      let rowhead = document.createElement("tr");
      Object.keys(pepTableColHeaders).forEach( n => { 
        let th = document.createElement('th'); 
        // th.className = 'peptb_'+n; 
        headname = pepTableColHeaders[n];
        if (headname=="Support") th.setAttribute("class","text-center");
        th.innerHTML = headname; 
        rowhead.appendChild(th); 
      });
      pepTablehead.appendChild(rowhead);
      pepTable.appendChild(pepTablehead);
      let pepTablebody = document.createElement("tbody");
      // let tr = document.createElement("tr");
      const thisrow = proteinrow["peptides"];
      for (let i=0; i<thisrow['Peptide'].length; i++) {
        let tr = document.createElement("tr");
        let tdpep = document.createElement("th"); tdpep.setAttribute("scope","row");
        let tdsupport = document.createElement("td");
        let tdt12 = document.createElement("td");
        let tdr2 = document.createElement("td");
        let thissupporti = thisrow['support'][i];
        let thist12i = thisrow['t12_CFit'][i];
        let thisr2i = (thissupporti>2)? thisrow['r2_CFit'][i]: "";

        tdpep.innerHTML = thisrow["Peptide"][i]; 
        tdsupport.innerHTML = thissupporti; tdsupport.setAttribute("class","text-center");
        tdt12.innerHTML = roundVal(thist12i,2);
        tdr2.innerHTML = roundVal(thisr2i,4);

        goodpep = (thissupporti>2 && thisr2i>0.8)?true:false;
        if (goodpep) { tr.setAttribute("class","table-warning"); } else { tr.style.fontStyle="italic"; } 
        tr.appendChild(tdpep);
        tr.appendChild(tdsupport); 
        tr.appendChild(tdt12);
        tr.appendChild(tdr2);
        pepTablebody.appendChild(tr);
      }
      pepTable.appendChild(pepTablebody);
      modalPeptideTb.replaceChild(pepTable, modalPeptideTb.firstChild) ; 
  
    })
  
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

function showProteinChart(elm, genename="") {
  const gene = (genename)?genename: elm.parentNode.id; // event.target.attributes['proteinlink'].value;

  // set focus 
  const position =  $("div#proteinChart iframe").position();
  scroll(0,position.top);

  var proteinframe = $("div#proteinChart iframe")[0];
  var peptideframe = $("div#peptideChart iframe")[0];
  const proteinPage = "./charts/proteinLevel/RelAbundance_Gene-"+gene+".html";
  const peptidePage = "./charts/peptideLevel/RelAbundance_Gene-"+gene+"-peptides.html";

  proteinframe.setAttribute('src', proteinPage);
  peptideframe.setAttribute('src', peptidePage);

  // change the pop-up links
  var proteinpopup = $("div#proteinChart a:first")[0];
  var peptidepopup = $("div#peptideChart a:first")[0];
  proteinpopup.setAttribute('href', proteinPage);
  peptidepopup.setAttribute('href', peptidePage);

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
