'use client';
import { Kekule } from 'kekule';
import { useEffect, useState } from "react";

export default function Board() {

  //Kekule.externalResourceManager.register('three.js', THREE);

  const [data, setData] = useState("")
  const [zoom, setZomm] = useState(2)
  const [currentTime, setCurrentTime] = useState(0);

  useEffect(() => {
    document.getElementById("name")?.focus()
  }, [])

  useEffect(() => {
    fetch('/api/time').then(res => res.json()).then(data => {
      setCurrentTime(data.time);
    });
  }, [])



  function draw() {

    const name = (document.getElementById("name") as HTMLInputElement).value;
    const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${name}/SDF`

    fetch(url).then((response) => {
      return response.ok ? response.text() : ""
    })
      .then((data) => {
        setData(data)
        if (data == "") {
          getDrawPanel()
        } else {
          var mol = readMolecule(data);
          var condensed = (document.getElementById('checkBoxCondensed') as HTMLInputElement).checked;
          drawMolecule(mol, condensed, zoom);
        }
      })
  }

  return (
    <div className="container mt-5">
      <h1 className="text-indigo-700 text-xl font-bold">BioActivitat 1 - Molecule Viewer</h1>
      <div className="grid grid-cols-2 gap-2">

        <div id="01">
          <input id="name" type="text" placeholder="adenine" onKeyUp={draw} />
          <div>
            <input className="form-check-input" type="checkbox" name="checkBoxCondensed" id="checkBoxCondensed" onClick={draw} />
            <label className="form-check-label ms-3" htmlFor="checkBoxCondensed">Condensed</label>
          </div>
          <div>
            <label>Zoom</label>
            <input type="number" min="1" max="5" value="2" />
          </div>
          <div id="drawStage border border-3 rounded">
            <div id="draw-panel"></div>
          </div>
        </div>

        <div id="02">
          <h6 className="text-center text-xl">SDF file</h6>
          <pre id="data-panel" className="p-3 border border-3 rounded text-wrap">
            {data}
          </pre>
        </div>
      </div>



      <p className="text-center mt-5">Kekule {Kekule.VERSION}</p>
      <p className="text-center mt-5">The current time is {currentTime}.</p>

    </div>
  );
}

var drawBoxWidth = 500;
var drawBoxHeight = 500;


function readMolecule(data: string) {

  var reader = new Kekule.IO.MdlReader();
  //var multiple = parseInt(document.getElementById('editMultiple').value);
  var mol1;
  var result;



  var r = reader.readData(data, Kekule.IO.ChemDataType.TEXT);
  mol1 = r;
  result = r;


  var ops = mol1.getRenderOptions();
  mol1.setRenderOption('atomColor', '#FF0000');
  mol1.setRenderOption('opacity', 1);
  if (mol1.hasCtab && mol1.hasCtab()) {
    var connector = mol1.getConnectorAt(0);
    if (connector !== undefined) {
      connector.setRenderOption('color', '#0000FF');
      connector.setRenderOption('opacity', 1);
    }
    var node = mol1.getNodeAt(0);
    node.setRenderOption('atomColor', '#000000');
    node.setRenderOption('opacity', 1);
  }

  return result;
}

function getDrawPanel(): HTMLElement {
  var drawPanel = document.getElementById('draw-panel') as HTMLElement;
  drawPanel.innerHTML = ""
  return drawPanel
}



function drawMolecule(mol: any, condensed: any, zoom: number) {

  var parentElem = getDrawPanel()

  var bridge = new Kekule.Render.CanvasRendererBridge();

  var context = bridge.createContext(parentElem, drawBoxWidth, drawBoxHeight);

  var baseCoord = { 'x': drawBoxWidth / 2, 'y': drawBoxHeight / 2 };
  var options = {
    moleculeDisplayType: condensed ? Kekule.Render.MoleculeDisplayType.CONDENSED : Kekule.Render.MoleculeDisplayType.SKELETAL,
    retainAspect: true,
    autoScale: true,
    zoom: zoom
  };

  bridge.clearContext(context);

  var painter = new Kekule.Render.ChemObjPainter(Kekule.CoordMode.COORD2D, mol, bridge);
  painter.draw(context, baseCoord, options);
}
