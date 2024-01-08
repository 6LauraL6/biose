import type { Metadata } from 'next'
import { Inter } from 'next/font/google'
import Navigation from './components/navigation/navbar'
import './globals.css'
import './locals.css'

const inter = Inter({ subsets: ['latin'] })

export const metadata: Metadata = {
  title: 'Bio - Sequence',
  description: '',
}

export default function RootLayout({
  children,
}: {
  children: React.ReactNode
}) {
  return (
    <html>
      <head>
        <title>Bio Sequence</title>
        <meta name="viewport" content="width=device-width, initial-scale=1" />
      </head>
      <body className={inter.className}>
        <Navigation />
        {children}
      </body>
    </html>
  )
}
