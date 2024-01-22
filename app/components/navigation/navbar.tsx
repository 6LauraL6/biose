// Navbar.tsx
'use client'// Navbar.tsx
import React, { useState } from "react";
import Link from "next/link";
import Logo from "./Logo";

const Navbar = () => {
  const [isOpen, setIsOpen] = useState(false);

  const toggleMenu = () => {
    setIsOpen(!isOpen);
  };

  const closeMenu = () => {
    setIsOpen(false);
  };

  return (
    <>
      <div className="w-full h-20 bg-indigo-300 sticky top-0">
        <div className="container px-4 h-full">
          <div className="flex justify-between items-center h-full">
            <div className="md:hidden">
              <button
                className="text-white focus:outline-none"
                onClick={toggleMenu}
              >
                <svg
                  className="w-6 h-6"
                  fill="none"
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth="2"
                  viewBox="0 0 24 24"
                  stroke="currentColor"
                >
                  <path d="M4 6h16M4 12h16m-7 6h7"></path>
                </svg>
              </button>
            </div>
            <ul
              className={`${
                isOpen ? "flex" : "hidden"
              } md:flex gap-x-6 text-white`}
            >
              <Logo />
              <li>
                <Link href="/act1" className="text-xl hover:text-yellow-500 transition hover:underline">
                  <p className="cursor-pointer" onClick={closeMenu} >
                    Pt1-Ses1-Kekule
                  </p>
                </Link>
              </li>
              <li>
                <Link href="/" className="text-xl hover:text-yellow-500 transition hover:underline">
                  <p className="cursor-pointer" onClick={closeMenu}>
                  Pt1-Ses2-DNA-RNA-Prot
                  </p>
                </Link>
              </li>
              <li>
                <Link href="/act4" className="text-xl hover:text-yellow-500 transition hover:underline">
                  <p className="cursor-pointer" onClick={closeMenu} >
                    Pt1-Ses3-Genbank
                  </p>
                </Link>
              </li>
            </ul>
          </div>
        </div>
      </div>
    </>
  );
};

export default Navbar;