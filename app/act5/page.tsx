'use client';
import { useEffect, useState } from "react";
import * as eutils from "ncbi-eutils";
import { Database } from 'react-bootstrap-icons';

export default function Home() {

  const [dbs, setDbs] = useState([])
  const [db, setDb] = useState(null)

  if (dbs.length == 0) {
    eutils.einfo().then(function (data: any) {
      setDbs(data.einforesult.dblist.sort())
    });
  }

  function dbInfo(name: string) {
    eutils.einfo(name).then(function (data: any) {
      setDb(data.einforesult.dbinfo[0])
    });
  }

  function home() {
    setDb(null)
  }


  var panel;
  if (db != null) {
    panel = <Db db={db} handle={home} />
  }
  else {
    if (dbs.length == 0) {
      const style = {
        width: "8rem", height: "8rem"
      };
      panel =
        <div className="text-center mt-5">
          <div className="spinner-border text-danger" style={style} role="status">
            <span className="sr-only">EInfo</span>
          </div>
        </div>;
    } else {
      panel = <>
        <p className="text-center">{dbs.length} databases found.</p>
        <div className="row mt-2">
          <div className="row mt-2">
            {dbs.map(item => (
              <div key={item} className="col-4 col-sm-2">
                <DbLink name={item} handle={dbInfo}></DbLink>
              </div>))
            }
          </div>
        </div>;
      </>
    };
  }

  return (
    <>
      <main className="container mt-5">
        <h1 className="text-center mt-5 fs-2">Entrez Databases</h1>
        <hr />
        {panel}
      </main>
    </>
  );
}


interface DbLinkProps {
  name: string
  handle: any
}
function DbLink({ name, handle }: DbLinkProps) {
  return <div className="p-2 m-2 border border-primary border-2 rounded text-center link-primary fw-bold" key={name} onClick={() => handle(name)} >{name}</div>
}


interface DbProps {
  db: any
  handle: any
}
function Db({ db, handle }: DbProps) {

  return <>
    <div className="row mt-2">
      <div className="col-2"><button type="button" onClick={() => handle()} className="btn btn-outline-primary">
        <Database/> Databases</button></div>
      <div className="col"><h2 className="fw-bold">{db.dbname}</h2></div>
    </div>
    <div className="border border-2 rounded-5 mt-2 border-success">
      <table className="table table-bordered ms-5 me-5 mt-3 mb-3">
        <tr><th className="font-monospace">menuname</th><td>{db.menuname}</td></tr>
        <tr><th className="font-monospace">description</th><td>{db.description}</td></tr>
        <tr><th className="font-monospace">count</th><td>{db.count}</td></tr>
        <tr><th className="font-monospace">lastupdate</th><td>{db.lastupdate}</td></tr>
        <tr><th className="font-monospace">dbbuild</th><td>{db.dbbuild}</td></tr>
      </table>
    </div>

    <h2 className="text-center mt-5">Linked List</h2>
    <table className="table table-bordered table-hover">
      <thead>
        <tr>
          <th>dbto</th>
          <th>name</th>
          <th>menu</th>
          <th>description</th>
        </tr>
      </thead>
      <tbody>
        {db.linklist.map((link: any) => (
          <tr key={link}><td>{link.dbto}</td><td>{link.name}</td><td>{link.menu}</td><td>{link.description}</td></tr>
        ))}
      </tbody>
    </table>


    <h2 className="text-center mt-5">Field List</h2>

    <table className="table table-bordered table-hover">
      <thead>
        <tr>
          <th>name</th>
          <th>fullname</th>
          <th>description</th>
          <th>termcount</th>
          <th>isdate</th>
          <th>isnumerical</th>
          <th>singletoken</th>
          <th>hierarchy</th>
          <th>ishidden</th>
          <th>istruncatable</th>
          <th>israngeable</th>
        </tr>
      </thead>
      <tbody>
        {db.fieldlist.map((field: any) => (
          <tr key={field}><td>{field.name}</td><td>{field.fullname}</td><td>{field.description}</td><td>{field.termcount}</td><td>{field.isdate}</td><td>{field.isnumerical}</td><td>{field.singletoken}</td><td>{field.hierarchy}</td><td>{field.ishidden}</td><td>{field.istruncable}</td><td>{field.rangeable}</td></tr>
        ))}
      </tbody>
    </table>
  </>
}